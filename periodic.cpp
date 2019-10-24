//#include "core.h"
//#include "matrix.h"
//#include "visualization.h"
//#include "model.h"
#include "lik.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "periodic.h"
#include "dynrec.h"
#include "ftsolver.h"

using namespace Eigen;

periodic::periodic(kinematicmodel* model_){
  model = model_;
  n_t = 0;
  traj_size = 0;
  traj = NULL;
  set_dynparts();
  dynrecs = NULL;
  ftsolver = new forcetorquesolver (this);
  joint_vel_traj = NULL;
  computed_torques = NULL;
}

periodic::~periodic(){
  delete_dynrecs();
  clear_traj();
  delete [] parentis;
  delete [] masses;
  delete [] footis;
  for(unsigned int i=0;i<dynparts.size();i++){delete dynparts[i];}
  delete ftsolver;
}

void periodic::set_dynparts(){
  map<modelnode*,int> mnode_id_map;
  set<modelnode*> foot_set;
  set_footset(foot_set);
  int id = 0, fi = 0;
  vector<odepart*>* odeparts = model->get_odeparts();
  int n = odeparts->size();
  parentis = new int [n];
  masses = new double [n];
  footis = new int [nfeet];
  vector<odepart*>::iterator it = odeparts->begin();
  for(;it!=odeparts->end();it++){
    odepart* opart = (*it);
    modelnode* mnode = opart->get_mnode();
    mnode_id_map[mnode] = id;
    dynpart* dpart = new dynpart (opart);
    dpart->setup(mnode_id_map);
    dpart->setup_foot(foot_set);
    if(dpart->if_foot()){footis[fi++] = id;}
    dynparts.push_back(dpart);
    parentis[id] = dpart->get_parent_id();
    masses[id] = dpart->get_mass();
    id++;
  }
}

void periodic::print(){
  cout << "--- periodic ---" << endl;
  cout << "number of dynparts: " << dynparts.size() << endl;
  cout << "dynparts:" << endl;
  vector<dynpart*>::iterator it = dynparts.begin();
  for(;it!=dynparts.end();it++){(*it)->print();}
  cout << "n_t = " << n_t << endl;
  cout << "config_dim = " << config_dim << endl;
  cout << "traj_size = " << traj_size << endl;
  cout << "dt_traj = " << dt_traj << endl;
  cout << "rcap = " << rcap << endl;
}

void periodic::record_trajectory(pergensetup* pgs, int n_t_){
  n_t = n_t_;
  traj_size = n_t + 5;
  config_dim = pgs->get_config_dim();
  nmj = model->number_of_motor_joints();
  double t = 0;
  //double dt = pgs->get_pergen()->get_period() / n_t;
  double dt = pgs->get_period() / n_t;
  double* rec = new double [config_dim];
  traj = new_2d_array(traj_size,config_dim);
  double** p = traj;
  for(int i=0;i<traj_size;i++){
    pgs->set_rec(rec,t);
    model->set_jvalues_with_lik(rec);
    model->get_jvalues(*p++);
    t += dt;
  }
  delete rec;
  dt_traj = dt;
  //rcap = pgs->get_rcap();
  rcap = model->get_lik()->get_rcap();
}

void periodic::print_trajectory(){
  cout << "--- trajectory (all) ---" << endl;
  cout << "traj_size = " << traj_size << endl;
  for(int i=0;i<traj_size;i++){
    print_trajectory(i);
  }
}

void periodic::print_trajectory(int i, int m){
  cout << "--- trajectory (" << i << " to " << (i+m-1)%traj_size <<") ---" << endl;
  cout << "traj_size = " << traj_size << endl;
  for(int j=0;j<m;j++){
    print_trajectory(i++);
  }
}

void periodic::print_trajectory(int i){
  i %= traj_size;
  double* p = traj[i];
  for(int j=0;j<config_dim;j++){
    if(j){cout << " ";}
    cout << *p++;
  }
  cout << endl;
}

void periodic::clear_traj(){
  delete_2d_array(traj,traj_size);
  delete_joint_vel_traj();
  delete_computed_torques();
  traj_size = 0;
}


void periodic::new_dynrecs(){
  dynrecs = new dynrecord* [traj_size];
  int size = dynparts.size();
  dynrecord** p_rec = dynrecs;
  for(int i=0;i<traj_size;i++){
    *p_rec++ = new dynrecord (size, nfeet);
  }
}

void periodic::delete_dynrecs(){
  dynrecord** p_rec = dynrecs;
  for(int i=0;i<traj_size;i++){delete *p_rec++;}
  delete [] dynrecs;
}

void periodic::compute_dynrecs(){
  if(!dynrecs){new_dynrecs();}
  double** p_traj = traj;
  dynrecord** p_rec = dynrecs;
  for(int i=0;i<traj_size;i++){
    model->set_jvalues(*p_traj++);
    model->recompute_modelnodes();
    //model->print(0);
    recompute_dynparts();
    (*p_rec++)->initialize(dynparts,rcap);
  }
}

void periodic::print_dynrecs(int id){
  cout << "--- dynrecs (all) ---" << endl;
  cout << "part id = " << id << endl; 
  for(int i=0;i<traj_size;i++){dynrecs[i]->print(id);}
}

void periodic::print_dynrecs(int id, int i, int m){
  cout << "--- dynrecs (" << i << " to " << (i+m-1)%traj_size <<") ---" << endl;
  cout << "part id = " << id << endl; 
  for(int j=0;j<m;j++){dynrecs[(i+j) % traj_size]->print(id);}
}

void periodic::print_dynrec(int i){
  cout << "--- dynrec (i = " << i <<") ---" << endl;
  for(int j=0;j<(int)dynparts.size();j++){
    cout << "part id = " << j << ":" << endl;
    dynrecs[i]->print(j);
  }
}

void periodic::recompute_dynparts(){
  vector<dynpart*>::iterator it = dynparts.begin();
  for(;it!=dynparts.end();it++){(*it)->recompute();}
}

// order n+1 ders are computed one stage after order n
void periodic::compute_dynrec_ders(){
  dynrecord** p = dynrecs;
  int n_stages = 2;
  for(int i=0;i<traj_size+n_stages;i++){
    for(int j=0;j<n_stages;j++){
      if(i > 2*j && i < traj_size-1){
	p[i-j]->compute_ders(j,p[i-1-j],p[i+1-j],dt_traj,this);
      }
    }
  }
}

void periodic::switch_torso_penalty(bool force, bool torque){
  ftsolver->switch_torso_penalty(force,torque);
}

void periodic::check_solve_ft(){
  //test_mats(); // temporary
  VectorXd x, y, xsum;
  double s=0;
  int i0=2;
  for(int i=i0;i<n_t+i0+1*0;i++){
    //solve_forcetorques(i,x,y);
    ftsolver->solve_forcetorques(dynrecs[i],x,y);
    if (i==i0) {xsum = x;}
    else {xsum += x;}
    //cout << "i="<<i<<" "<<x.transpose().head(3) << endl;
    //cout<<y.transpose()<<endl;
    //for(int j=2;j<y.size();j+=3){s += y(j);}
    cout<<"z: ";for(int j=2;j<y.size();j+=3){cout<<y(j)<<" ";}cout<<endl;
    s+=y.sum();
  }
  x = xsum / n_t;
  s /= n_t;
  VectorXd x1 = x.head(y.size()/2);
  //cout << "l1 ynorm="<<l1_norm(y1) << endl;
  cout << "l1 ynorm = " << x1.lpNorm<1>() << endl;
  cout << "xnorm="<<x.head(x.size()/2).norm() << endl;
  cout << "xnorm="<<x.norm() << endl;
  //cout << y.head(10) << endl;
  cout << x.head(x.size()/2) << endl;
  cout<<"s = "<<s<<endl;
}

void periodic::set_footset(set<modelnode*>& foot_set){
  model->get_foot_mnodes(foot_set);
  nfeet = foot_set.size();
}

void periodic::new_joint_vel_traj(){
  joint_vel_traj = new_2d_array(traj_size,config_dim);
}

void periodic::delete_joint_vel_traj(){
  delete_2d_array(joint_vel_traj,traj_size);
}

void periodic::new_computed_torques(){
  //nmj = model->number_of_motor_joints();
  computed_torques = new_2d_array(n_t,nmj);
}

void periodic::delete_computed_torques(){
  delete_2d_array(computed_torques,n_t);
}

void periodic::compute_joint_vel_traj(){
  double** p = joint_vel_traj;
  if(p == NULL){
    new_joint_vel_traj();
    p = joint_vel_traj;
  }
  p++;
  double** p1 = traj;
  double** p2 = p1+2;
  double *p_, *p1_, *p2_;
  for(int i=2;i<traj_size;i++){
    p_ = *p++;
    p1_ = *p1++;
    p2_ = *p2++;
    for(int j=0;j<config_dim;j++){
      double d = (*p2_++ - *p1_++);
      if (d > M_PI) {d -= 2*M_PI;}
      else if (d < -M_PI) {d += 2*M_PI;}
      *p_++ = d/(2*dt_traj);
    }
  }
}

double periodic::work_over_period(){
  compute_torques_over_period();
  compute_joint_vel_traj();
  //nmj = model->number_of_motor_joints();
  int i0 = 2;
  double** p_jvtraj = joint_vel_traj + i0;
  double work_period = 0; // work over period
  for(int i=i0;i<n_t+i0;i++){
    double work_dt = 0; // work over dt
    double *p = (*p_jvtraj++) + 6, *p1 = computed_torques[i % n_t];
    for(int j=0;j<nmj;j++){
      double torqa = *p1++;
      double joint_vel = *p++;
      double dw = torqa*joint_vel;
      //dw = fabs(dw); // counts |dw|, not just positive dw
      dw = (dw > 0)? dw : 0;
      work_dt += dw;
    }
    work_dt *= dt_traj;
    work_period += work_dt;
  }
  //cout << "work over period = " << work_period << endl;
  return work_period;
}

void periodic::hinge_joint_part_ids(list<int>& ids){
  int n = get_number_of_dynparts();
  for(int i=0;i<n;i++){
    modeljoint* joint = model->get_mnode(i)->get_joint();
    if(joint){
      if(joint->get_type() == hinge){ids.push_back(i);}
    }
  }
}

double periodic::get_total_mass(){
  int n = get_number_of_dynparts();
  double m_tot = 0;
  for(int i=0;i<n;i++){m_tot += masses[i];}
  return m_tot;
}

void periodic::get_motor_torques(double* motor_torques){
  int n = get_number_of_dynparts();
  list<int> ids; // ids of parts with joints
  hinge_joint_part_ids(ids);
  double* p = motor_torques;
  VectorXd x = *ftsolver->get_fts();
  double* jzaxis = ftsolver->get_jzaxis();
  list<int>::iterator it = ids.begin();
  for(;it!=ids.end();it++){
    int k = 3*(*it);
    int k1 = 3*n+k;
    double motor_torque = 0;
    for(int j=0;j<3;j++){motor_torque += jzaxis[k+j]*x(k1+j);}
    *p++ = motor_torque;
  }
}

void periodic::analyze_contforces(double* contforces){
  double *p = contforces;
  for(int i=0;i<nfeet;i++){
    double cfx, cfy, cfz;
    cfx = *p++; cfy = *p++; cfz = *p++;
    if(cfz < min_cfz){min_cfz = cfz;}
    double cfxy = sqrt(cfx*cfx+cfy*cfy);
    double mu = cfxy/cfz;
    if(mu > max_mu){max_mu = mu;}
  }
}

// solves for motor torques and contact forces 
// given trajectory and contacts info
void periodic::solve_torques_contforces(int i, double* torques, double* contforces){
  VectorXd x, y;
  ftsolver->solve_forcetorques(dynrecs[i],x,y);
  get_motor_torques(torques);
  VectorXd::Map(contforces,y.size()) = y;
}

// solves for contact forces given trajectory and motor torques
void periodic::solve_contforces_given_torques(int i, double* contforces, double* torques){
  VectorXd y;
  //nmj = model->number_of_motor_joints();
  VectorXd z = Map<VectorXd> (torques,nmj);
  ftsolver->solve_forces(dynrecs[i],z,y);
  VectorXd::Map(contforces,y.size()) = y;
}

void periodic::compute_torques_over_period(){
  new_computed_torques();
  double* contforces = new double [3*nfeet];
  min_cfz = 1e10, max_mu = -1e10;
  int i0 = 2;
  for(int i=i0;i<n_t+i0;i++){
    //dynrecord* dynrec = dynrecs[i];
    VectorXd x, y;
    ftsolver->solve_forcetorques(dynrecs[i],x,y);
    VectorXd::Map(contforces,y.size()) = y;
    analyze_contforces(contforces);
    double* motor_torques = computed_torques[i % n_t];
    get_motor_torques(motor_torques);
  }
  delete [] contforces;
}

void periodic::get_motor_adas(int tsi, double* as, double* das){
  if(!(traj && joint_vel_traj)){cout<<"ERROR: no data"<<endl;exit(1);}
  tsi %= n_t;
  if (tsi < 2) {tsi += n_t;}
  int offset = 6;
  double *p = traj[tsi]+offset, *p1 = joint_vel_traj[tsi]+offset;
  for(int i=offset;i<config_dim;i++){
    *as++ = *p++;
    *das++ = *p1++;
  }
}

void periodic::get_complete_traj_rec(int tsi, double* rec){
  if(tsi >= n_t){cout<<"ERROR: time step must be < n_t"<<endl;exit(1);}
  if (tsi < 2) {tsi += n_t;}
  arrayops ao (config_dim), ao1 (nmj);
  ao.assign(rec,traj[tsi]);
  rec += config_dim;
  ao.assign(rec,joint_vel_traj[tsi]);
  rec += config_dim;
  ao1.assign(rec,get_computed_torques(tsi));
}

void periodic::get_complete_traj(double** complete_traj){
  compute_torques_over_period();
  compute_joint_vel_traj();
  double **p = complete_traj;
  for(int i=0;i<n_t;i++){
    get_complete_traj_rec(i,*p++);
  }
}
