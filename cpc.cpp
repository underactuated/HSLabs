//#include "core.h"
//#include "matrix.h"
//#include "visualization.h"
//#include "model.h"
//#include "pergen.h"
//#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include "periodic.h"
//#include "dynrec.h"
#include "cpc.h"
#include "effdata.h"

cand::cand(int tpi_, double t0_, double s_, double loss_){
  tpi = tpi_; 
  t0 = t0_; 
  s = s_; 
  loss = loss_;
}

void cand::print(){
  cout << "tpi = " << tpi << ", t0 = " << t0 << ", s = " << s << endl;
}


cpccontroller::cpccontroller(const kinematicmodel* model){
  q_dim = model->get_config_dim();
  chi_dim = model->number_of_motor_joints();
  psi_dim = q_dim - chi_dim;
  B.resize(q_dim,chi_dim);
  Bt.resize(chi_dim,q_dim);
  b.resize(q_dim,psi_dim);
  target_points = NULL;
  tps_size = 0;
  current_state = new double [2*q_dim];
  last_cand.tpi = 0;
  effdata = NULL;
  set_flag("poscontrol",false);
  set_flag("goal_t0s",false);
  set_flag("tpdist_switch",false); // swithes between cpc s and s=sg
  int maska[] = {0,1,5};
  mask.insert(mask.end(), maska, maska+3);
  low_tpdist = true;
}

cpccontroller::~cpccontroller(){
  if(target_points){delete_2d_array(target_points,tps_size);}
  delete [] current_state;
  if(effdata){delete effdata;}
}

void cpccontroller::set_target_points_by_per(periodic* per){
  //per->print();exit(1);
  if(target_points){cout<<"ERROR: target points are present"<<endl;exit(1);}
  int nt = per->get_nt();
  int rec_len = 2*q_dim + chi_dim;
  target_points = new_2d_array(nt,rec_len);
  for(int i=0;i<nt;i++){
    per->get_complete_traj_rec(i,target_points[i]);
    apply_mask(target_points[i]);
  }
  tps_size = nt;
  //print_target_points();exit(1);
}

// we probably do not need B, only Bt
void cpccontroller::set_B_by_player(double** B_transp){
  for(int i=0;i<chi_dim;i++){
    B.col(i) = Map<VectorXd> (B_transp[i],q_dim);
    Bt.row(i) = Map<VectorXd> (B_transp[i],q_dim);
  }
  compute_b();
}

void cpccontroller::print(){print(0);}

void cpccontroller::print(int detail_level){
  cout << "----- CPC controller -----" << endl;
  cout << "w = " << w << ", sg = " << sg;
  cout << ", k0 = " << k0 << ", kc = " << kc;
  cout << ", tauc = " << tauc << ", n_cand = " << n_cand << endl;
  cout << "q_dim = " << q_dim << ", chi_dim = " << chi_dim;
  cout << ", psi_dim = " << psi_dim << endl;
  cout << "tps_size = " << tps_size << endl;

  if(detail_level>0){
    /*cout << "target points:" << endl;
    int rec_len = 2*q_dim+chi_dim;
    for(int i=0;i<tps_size;i++){
      print_array<double>(target_points[i],rec_len);
      }*/
    print_target_points();
  }
}

void cpccontroller::print_target_points(){
  cout << "----- target points -----" << endl;
  int rec_len = 2*q_dim + chi_dim;
  for(int i=0;i<tps_size;i++){
    print_array(target_points[i],rec_len);
  }
}

void cpccontroller::set_current_state(const double* config, const double* dconfig){
  arrayops ao (q_dim);
  ao.assign(current_state,config);
  ao.assign(current_state+q_dim,dconfig);
  apply_mask(current_state);
}

//int tmp=0;
void cpccontroller::compute_b(){
  MatrixXd Bt_chi = Bt.rightCols(chi_dim);
  B_chi = Bt_chi.transpose();
  b.topRows(psi_dim) = -MatrixXd::Identity(psi_dim,psi_dim);

  ColPivHouseholderQR<MatrixXd> decomp(Bt_chi);
  for(int i=0;i<psi_dim;i++){
    b.block(psi_dim,0,chi_dim,psi_dim).col(i) = decomp.solve(Bt.col(i));
  }

  bbt = b*b.transpose();

  /*b0 = b;
  b0.topRows(2) *= 0;
  b0b0t = b0*b0.transpose();*/

  //if(rand()%10==0){cout<<b<<endl;exit(1);}
  //cout<<b<<endl;exit(1);

  //compute_prox_loss_over_tpset();
}

void cpccontroller::compute_prox_loss_over_tpset(){
  list<int> tpis;
  map<double,cand> loss_cands;
  for(int i=0;i<tps_size;i++){tpis.push_back(i);}
  compute_prox_loss_over_tps(tpis,loss_cands);
  //exit(1);
}

void cpccontroller::state_to_qdq(double* state, VectorXd& q, VectorXd& dq){
  q = Map<VectorXd> (state,q_dim);
  dq = Map<VectorXd> (state+q_dim,q_dim);
}

void cpccontroller::state_to_chidchi(double* state, VectorXd& chi, VectorXd& dchi){
  chi = Map<VectorXd> (state+psi_dim,chi_dim);
  dchi = Map<VectorXd> (state+q_dim+psi_dim,chi_dim);
}

void cpccontroller::compute_prox_loss_over_tps(const list<int>& tpis, map<double,cand>& loss_cands){
  VectorXd qd, dqd, q0, dq0;
  state_to_qdq(current_state,q0,dq0);
  MatrixXd bt = b.transpose();
  //VectorXd btildt = (bbt*dq0).transpose();// experim

  list<int>::const_iterator it = tpis.begin();
  for(;it!=tpis.end();it++){
    int tpi = (*it);
    state_to_qdq(target_points[tpi],qd,dqd);
    VectorXd btildt = (bbt*dqd).transpose();

    //qd.segment(0,2)*=0;dqd.segment(0,2)*=0;q0.segment(0,2)*=0;dq0.segment(0,2)*=0; // experim
    
    VectorXd del_q = qd - q0;
    //del_q.segment(0,2) *= 0;
    //del_q.segment(0,1) *= 0;
    double denom = btildt.dot(dq0);
    double t0 = btildt.dot(del_q)/denom;
    double s = btildt.dot(dqd)/denom;

    VectorXd c0 = bt*(del_q-dqd*t0/s);
    VectorXd c1 = bt*(dqd/s-dq0);
    double cc0 = c0.dot(c0), cc1 = c1.dot(c1);

    double loss = prox_loss(t0,s,cc0,cc1);
    loss_cands[loss] = cand (tpi,t0,s,loss);

    //cout << "t0 = " << t0 << ", s = " << s << ", cc0 = " << cc0 << ", cc1 = " << cc1 << ", loss = " << loss << endl;
  }
}

double cpccontroller::prox_loss(double t0, double s, double cc0, double cc1){
  double a = w*t0, a1 = s-sg;
  //double a = w*t0, a1 = 1./s-1./sg; // experim
  //if(s<0.5){a1*=10000;}
  //if(fabs(a1)>.5){a1*=100000;}
  return a*a + a1*a1+10*(cc0+cc1); // 100*cc for kicks
  //double ssg = s*sg; return a*a + a1*a1*(1./(ssg*ssg)+1.)+100*(cc0+cc1);
}

void cpccontroller::set_w_sg_nd_k0_kc_tauc(double w_, double sg_, int nd, double k0_, double kc_, double tauc_){
  w = w_;
  sg = sg_;
  n_cand = nd;
  k0 = k0_;
  kc = kc_;
  tauc = tauc_;
}

void clip_norm(VectorXd& tau, double tauc){
  double norm = tau.norm();
  if(norm > tauc){tau *= (tauc/norm);}
}

double smallest_eigenvalue(MatrixXd& m){
  VectorXd evs = m.eigenvalues().transpose().cwiseAbs();
  double ev_min = 1e10;
  for(int i=0;i<evs.size();i++){
    if(evs(i)<ev_min){ev_min = evs(i);}
  }
  return ev_min;
}

void cpccontroller::get_motor_torques(double* torques){
  //cout<<B_chi.eigenvalues().transpose()<<endl;
  //cout<<smallest_eigenvalue(B_chi)<<endl;
  //cout<<B_chi.determinant()<<endl;//cin.ignore();
  //print_array<double>(current_state+6,6);print_array<double>(target_points[last_cand.tpi]+6,6);cin.ignore();
  if(poscontrol_flag){//cout<<B_chi.eigenvalues().transpose()<<endl;
    B_chi_dec.compute(B_chi);
    get_poscontrol_torques(torques);
    return;
  }
  list<cand> cands;
  candidates(cands);
  B_chi_dec.compute(B_chi);
  MatrixXd cand_taus (chi_dim,n_cand);
  double k = k0;
  int i_best;
  do {
    //i_best = cost(k,cands,cand_taus);
    i_best = cost2(k,cands,cand_taus);
    //i_best = cost1(k,cands,cand_taus);
  //controls(torques,cands.front());
    k /= 2;
  } while (k > kc && cand_taus.col(i_best).norm() > tauc);
  //cout<<"k/k0 = "<<2*k/k0<<endl;
  vector<cand> candsv (cands.begin(),cands.end());
  last_cand = candsv[i_best];
  //cout<<last_cand.s<<endl;
  //last_cand.print();
  tp_dist_check();

  VectorXd tau = cand_taus.col(i_best);
  clip_norm(tau,tauc);
  VectorXd::Map(torques,chi_dim) = tau;

  //cout<<cands.size()<<endl;exit(1);
}

void cpccontroller::candidates(list<cand>& cands){
  list<int> tpis;
  map<double,cand> loss_cands;
  tpis_subset(tpis);
  compute_prox_loss_over_tps(tpis,loss_cands);
  map<double,cand>::iterator it = loss_cands.begin();
  for(int i=0;i<n_cand;i++){
    //cout << (*it).first << " " << (*it).second.tpi << endl;
    //double s = (*it).second.s; s=(s<1)?1:s; (*it).second.s=s; // experim
    if(goal_t0s_flag){
      cand* can = &((*it).second); 
      //can->t0=0; // setting only s=sg may be slightly better
      if (tpdist_switch_flag) {
	if (low_tpdist) {can->s = sg;}
      } else {can->s = sg;}
      //can->s=sg;
    }
    cands.push_back((*it).second);
    it++;
  }
}

void cpccontroller::controls(double* torques, double k, const cand& can){
  VectorXd tau;
  controls(tau,k,can);
  VectorXd::Map(torques,tau.size()) = tau;
}

void cpccontroller::controls(VectorXd& tau, double k, const cand& can){
  int tpi = can.tpi;
  double t0 = can.t0, s = can.s;

  VectorXd chi0, dchi0, chid, dchid;
  state_to_chidchi(current_state,chi0,dchi0);
  state_to_chidchi(target_points[tpi],chid,dchid);

  VectorXd Kdelchi = k*(chi0-chid+dchid*t0/s) + 2*sqrt(k)*(dchi0-dchid/s);
  VectorXd del_tau = -B_chi_dec.solve(Kdelchi);
  get_tau(tpi,tau);
  tau += del_tau;
  //VectorXd::Map(torques,tau.size()) = tau;
}

// same as controls, but has del_tau arg
void cpccontroller::controls1(VectorXd& tau, VectorXd& del_tau, double k, const cand& can){
  int tpi = can.tpi;
  double t0 = can.t0, s = can.s;

  VectorXd chi0, dchi0, chid, dchid;
  state_to_chidchi(current_state,chi0,dchi0);
  state_to_chidchi(target_points[tpi],chid,dchid);

  VectorXd Kdelchi = k*(chi0-chid+dchid*t0/s) + 2*sqrt(k)*(dchi0-dchid/s);
  del_tau = -B_chi_dec.solve(Kdelchi);
  get_tau(tpi,tau);
  tau += del_tau;
}

void cpccontroller::get_tau(int tpi, VectorXd& tau){
  tau = Map<VectorXd> (target_points[tpi]+2*q_dim,chi_dim);
}

int cpccontroller::cost(double k, const list<cand>& cands, MatrixXd& cand_taus){
  //cout<<"cost ---------------------"<<endl;
  //cout<<"k="<<k<<endl;
  list<cand>::const_iterator it = cands.begin();
  double norm_min = 1e10;
  int i=0, i_min=-1;
  for(;it!=cands.end();it++){
    VectorXd tau;
    controls(tau,k,*it);
    cand_taus.col(i) = tau;
    double norm = tau.norm();
    if(norm < norm_min){i_min = i; norm_min = norm;}
    //cout << norm << endl;
    i++;
  }
  return i_min;
}

double vec_mean(const vector<double>& vec){
  double s = 0;
  vector<double>::const_iterator it = vec.begin();
  for(;it!=vec.end();it++){s += (*it);}
  return s/vec.size();
}

int cpccontroller::cost2(double k, const list<cand>& cands, MatrixXd& cand_taus){
  vector<double> norms (n_cand), losss (n_cand);
  list<cand>::const_iterator it = cands.begin();
  int i = 0;
  for(;it!=cands.end();it++){
    losss[i] = (*it).loss;
    VectorXd tau;
    controls(tau,k,*it);
    cand_taus.col(i) = tau;
    norms[i] = tau.norm();
    i++;
  }

  double mean_norm = vec_mean(norms);
  double mean_loss = vec_mean(losss);

  double score_min = 1e10;
  int i_min=-1;
  for(int i=0;i<n_cand;i++){
    double score = norms[i]/mean_norm + losss[i]/mean_loss;
    if(score < score_min){i_min = i; score_min = score;}
  }
  return i_min;
}

// uses del_tau
int cpccontroller::cost1(double k, const list<cand>& cands, MatrixXd& cand_taus){
  //cout<<"cost ---------------------"<<endl;
  list<cand>::const_iterator it = cands.begin();
  double norm_min = 1e10;
  int i=0, i_min=-1;
  for(;it!=cands.end();it++){
    VectorXd tau, del_tau;
    controls1(tau,del_tau,k,*it);
    cand_taus.col(i) = tau;
    double norm = del_tau.norm();
    if(norm < norm_min){i_min = i; norm_min = norm;}
    //cout << tau.norm() << endl;
    i++;
  }
  return i_min;
}

void cpccontroller::load_tpset(string fname){
  if(target_points){cout<<"ERROR: target points are present"<<endl;exit(1);}
  ifstream file;
  file.open(fname.c_str());
  int rec_len = 2*q_dim + chi_dim;
  vector<double*> recs;
  string str;
  while(getline(file,str)){
    double* rec = new double [rec_len];
    stringstream ss; ss << str;
    //for(int i=0;i<rec_len;i++){ss >> rec[i];}
    for(int i=0;i<rec_len;i++){
      if(!(ss >> rec[i])){
	cout<<"ERROR: str too short in load_tpset"<<endl;exit(1);
      }
    }
    recs.push_back(rec);
    apply_mask(rec);
  }
  tps_size = recs.size();
  target_points = new_2d_array(tps_size,rec_len);
  arrayops ao (rec_len);
  for(int i=0;i<tps_size;i++){
    ao.assign(target_points[i],recs[i]);
    delete [] recs[i];
  }
  recs.clear();
  //print_target_points();exit(1);
}

void cpccontroller::get_torques0(double* torques){
  VectorXd tau;
  get_tau(last_cand.tpi,tau);
  //tau *= 0; cout<<tau<<endl;exit(1);
  VectorXd::Map(torques,tau.size()) = tau;
}

void cpccontroller::get_poscontrol_torques(double* torques){
  VectorXd tau;
  double k = k0;
  cand can (poscontrol_tpi++,0,1,-1);
  /*controls(tau,k,can);
  double norm = tau.norm();
  double tauc1 = 1*tauc;
  if(norm > tauc1){tau *= (tauc1/norm);}*/

  k = k0;
  do {controls(tau,k,can); k /= 2;
  } while (k > kc && tau.norm() > tauc);
  clip_norm(tau,tauc);

  VectorXd::Map(torques,tau.size()) = tau;
  poscontrol_tpi %= tps_size;
}

void cpccontroller::tpis_subset(list<int>& tpis){
  if(effdata){
    tpis_effdata(tpis,n_cand+100);
    return;
  }
  for(int i=0;i<tps_size;i++){tpis.push_back(i);}
}

void cpccontroller::set_flag(string flag_name, bool value){
  if(flag_name == "poscontrol"){
    poscontrol_flag = value;
    poscontrol_tpi = 0;
  } else if(flag_name == "goal_t0s"){
      goal_t0s_flag = value;
      set_flag("tpdist_switch",true);
  } else if(flag_name == "tpdist_switch"){
      tpdist_switch_flag = value;
  } else {
    cout << "ERROR: unknown flag " << flag_name << " in cpccontroller::set_flag" << endl; exit(1);
  }
}

void cpccontroller::setup_effdata(){
  if(effdata == NULL){effdata = new efficientdata;}
  effdata->prepare_data(target_points, tps_size, q_dim);
}

void cpccontroller::tpis_effdata(list<int>& tpis, int n){
  VectorXd q0, dq0;
  state_to_qdq(current_state,q0,dq0);
  //q0.segment(0,2) *= 0;
  VectorXd btil0 = bbt*dq0;
  double btq0 = btil0.dot(q0), btdq0 = btil0.dot(dq0); 
  double* btil = new double [q_dim];
  VectorXd::Map(btil, q_dim) = btil0;
  effdata->get_tpis(tpis, n, btil, btq0, btdq0, sg, w);
  delete [] btil;
}

void cpccontroller::apply_mask(double* state){
  vector<int>::const_iterator it = mask.begin();
  for(;it!=mask.end();it++){
    int i = (*it);
    state[i] = 0;
    state[q_dim + i] = 0;
  }
}

void cpccontroller::tp_dist_check(){
  int tpi = last_cand.tpi;
  arrayops ao (2*q_dim);
  //arrayops ao (q_dim);
  double dist = ao.distance(current_state,target_points[tpi]);
  low_tpdist = (dist < 0.5);
  //low_tpdist = (dist < 0.2);
  //cout<<"dist="<<dist<<endl;
  //exit(1);
} 

// Fits a plane to points my minimizing square error of z components.
// Arg plane is set to (c0,c1,c2), so it defines the fitted plane by
// defining its z component as z = c0*x+c1*y+c2.
void fit_plane(extvec& plane, list<extvec>& points){
  int n = points.size();
  //cout<<"n="<<n<<endl;
  MatrixXd m (n,3);
  VectorXd b (n);
  list<extvec>::iterator it = points.begin();
  for(int i=0;i<n;i++){
    extvec* point = &(*it++);
    double* p = point->get_data();
    m(i,0) = *p++;
    m(i,1) = *p++;
    m(i,2) = 1;
    b(i) = *p;
  }
  VectorXd x = m.colPivHouseholderQr().solve(b);
  /*cout << m << endl;
  cout << b.transpose() << endl;
  cout << x.transpose() << endl;*/
  plane.set(x(0),x(1),x(2));
  //plane.print();exit(1);
}

// TODO: old version, to remove
void fit_plane1(extvec& plane, list<extvec>& points){
  //return fit_plane1(plane,points);
  int n = points.size();
  //cout<<"n="<<n<<endl;
  MatrixXd m (n,3);
  list<extvec>::iterator it = points.begin();
  for(int i=0;i<n;i++){
    extvec* point = &(*it++);
    m.row(i) = Map<VectorXd> (point->get_data(),3);
  }
  VectorXd b = VectorXd::Constant(n,1);
  VectorXd x = m.colPivHouseholderQr().solve(b);
  /*cout << m << endl;
  cout << b.transpose() << endl;
  cout << x.transpose() << endl;*/
  plane.set(x(0),x(1),x(2));
  //plane.print();exit(1);
}
