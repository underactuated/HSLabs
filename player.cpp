#include "lik.h"
#include "pergen.h"
#include "player.h"
#include "periodic.h"
#include "odestate.h"
#include "cpc.h"
#include "geom.h"
#include "ghost.h"

// Step function is called from visualizer loop.
void visualizer::step(){
  if(player){player->step();}
}


modelplayer::modelplayer(){
  model = new kinematicmodel (true);
  visualizer* vis = get_vis(); // vis is initialized in visualization.cpp
  vis->set_player(this);
  set_default_flags();
  play_t = 0;
  play_dt = .01; // default dt
  play_rec = NULL;
  last_motor_torques = NULL;
  tdertrack = NULL;
  B_transp = NULL;
  hfield = NULL;
  ghost = NULL;
}

modelplayer::~modelplayer(){
  delete model;
  if(play_rec){delete [] play_rec;}
  if(last_motor_torques){delete [] last_motor_torques;}
  if(tdertrack){delete tdertrack;}
  if(B_transp){delete_2d_array(B_transp,nmj);}
  if(hfield){delete hfield;}
  if(ghost){delete ghost;}
}

void modelplayer::load_model(string fname){
  model->load_fromxml(fname);
  model->orient_odebodys(); 
  // ode-bodies must be oriented before setting ode-joints.
  model->set_ode_joints();
  make_play_rec();
}

bool cin_flag=false;//true;
// Step is called repeatedly from vis loop.
void modelplayer::step(){//cout<<"t = "<<play_t<<endl;
  if(cin_flag){cin.ignore();cin_flag=false;}
  //if(play_t > .9){cin.ignore();}
  switch (step_mode) {
  case 0: test0(); break; // myant.xml test
  case 1: test1(); break; // myant.xml test
  case 2: test2(); break; // hexapod.xml test
  case 3: test3(); break; // myant.xml test
  case 4: test4(); break; // myant.xml test
  case 5: play_pergensu_step(); break;
  case 6: simulate_ode(); return; break;
  case 7: break;
  default: 
    cout << "ERROR: case " << step_mode << " undefined in step" << endl; exit(1);
    break;
  }
  model->recompute_modelnodes();
  model->orient_odebodys();
}

// various (mostly old) tests, testi in [0,4]
void modelplayer::test(int testi){
  step_mode = testi;
  model->print();
  model->draw();
}

void modelplayer::set_flag(string flag_name, bool value){
  if(flag_name == "manual_viewpoint"){
    manual_viewpoint_flag = value;
    get_vis()->set_flag("manual_viewpoint",value);
  } else if(flag_name == "contact_force"){
    contact_force_flag = value;
  } else if(flag_name == "open_loop"){
    open_loop_flag = value;
  } else if(flag_name == "position_control"){
    position_control_flag = value;
  } else if(flag_name == "texture"){
    get_vis()->set_flag("texture",value);
  } else if(flag_name == "smooth_view"){
    get_vis()->set_flag("smooth_view",value);
  } else if(flag_name == "dynamics_from_simulation"){
    dynamics_from_simulation_flag = value;
    if(value){
      tdertrack = new configtimedertrack (2,config_dim,play_dt);
      B_transp = new_2d_array(nmj,config_dim);
    }
  } else if(flag_name == "cpc_control"){
    cpc_control_flag = value;
  } else if(flag_name == "record_traj"){
    record_traj_flag = value;
  } else if(flag_name == "torso_kicks"){
    torso_kicks_flag = value;
  } else if(flag_name == "ghost_walking"){
    ghost_walking_flag = value;
  } else if(flag_name == "clip_torque"){
    clip_torque_flag = value;
  } else {
    cout << "ERROR: unknown flag " << flag_name << " in modelplayer::set_flag" << endl; exit(1);
  }
}

// set_jangles_directly. do we need it?

// Sets joint angles on the model using LIK.
// rec contains torso orientation and foot positions.
void modelplayer::set_jangles_with_lik(const double* rec){
  model->set_jvalues_with_lik(rec);
}

void modelplayer::print_hip_poss(){
  model->get_lik()->print_hip_poss();
}

// partial setup of a pgs (pgs = pergen setup, pergen = periodic generator)
// Sets up pgs for given orientation and duration. Must be
// followed with setting TLh (period, step length and height).
void modelplayer::partial_setup_pergen(pergensetup& pergensu, const extvec* orientation, double step_duration){
  check_model_loaded();
  pgssweeper sweeper (NULL,model);
  sweeper.partial_setup_pergen(pergensu, orientation, step_duration);
}

// full setup of a pgs, for given pgsconfigparams
void modelplayer::setup_pergen(pergensetup& pergensu, const pgsconfigparams& pcp){
  check_model_loaded();
  pgssweeper sweeper (NULL,model);
  sweeper.setup_pergen(pergensu, &pcp);
}

// Makes a pergen setup according to setup specifications
// in a config_fname file.
pergensetup* modelplayer::make_pergensu(string config_fname, int setup_id){
  string confrec;
  get_rec_str(confrec,config_fname,setup_id);
  //cout << confrec << endl;
  pgsconfigparams pcp;
  get_pgs_config_params(confrec,pcp);
  string model_fname;
  model_fname = model->get_xmlfname();
  if (model_fname == "") {
    load_model(pcp.fname);
  } else {
    if (model_fname != pcp.fname) {
      cout<<"ERROR: model not from "<<pcp.fname<<endl;exit(1);
    }
  }
  int n = model->get_lik()->get_number_of_limbs();
  pergensetup* pergensu = new pergensetup (n);
  setup_pergen(*pergensu, pcp);
  return pergensu;
}

// Reads parameters from a pgs configuration string.
// Stores them in pgsconfigparams structure.
void modelplayer::get_pgs_config_params(const string& rec_str, pgsconfigparams& pcp){
  stringstream ss;
  ss << rec_str;
  extvec torso_pos, torso_angles;
  double period, step_length, step_height;
  string key;
  while(ss >> key){
    if (key == "xml_file") {ss >> pcp.fname;} 
    else if (key == "torso_pos") {
      double x, y, z;
      ss >> x >> y >> z;
      torso_pos.set(x,y,z);
    } else if (key == "torso_angles") {
      double x, y, z;
      ss >> x >> y >> z;
      torso_angles.set(x,y,z);
    } else if (key == "step_duration") {ss >> pcp.step_duration;} 
    else if (key == "period") {ss >> period;} 
    else if (key == "step_length") {ss >> step_length;}
    else if (key == "step_height") {ss >> step_height;} 
    else if (key == "curvature") {ss >> pcp.curvature;} 
    else if (key == "lateral_foot_shift") {
      double shift;
      ss >> shift;
      pcp.foot_shift = pair<int,double> (0,shift);
    } else if (key == "radial_foot_shift") {
      double shift;
      ss >> shift;
      pcp.foot_shift = pair<int,double> (1,shift);
    } else {
      cout<<"ERROR: unknown key "<<key<<endl;exit(1);
    }
  }
  pcp.set_TLh(period, step_length, step_height);
  //pcp.orientation[0].copy(torso_pos);
  //pcp.orientation[1].copy(torso_angles);
  pcp.orientation[0] = torso_pos;
  pcp.orientation[1] = torso_angles;
}

// Plays pergen setup from play_t = 0.
// This is not a simulation, the played trajectory
// need not be physically realizable.
void modelplayer::play_pergensu(pergensetup* pgs){
  play_pgs = pgs;
  play_t = 0;
  //set_flag("manual_viewpoint",false);
  step_mode = 5;
  model->draw();
}

// one step updates to play play_pgs
void modelplayer::play_pergensu_step(){
  play_pgs->set_rec(play_rec,play_t);
  set_jangles_with_lik(play_rec);
  play_t += play_dt;
}

// Gets a record, prefexed with "id ", from file fname,
// where id is int.
void modelplayer::get_rec_str(string& rec_str, string fname, int rec_id){
  ifstream file;
  file.open(fname.c_str());
  string str;
  while(getline(file,str)){
    stringstream ss; ss << str;
    int id;
    ss >> id;
    if(id == rec_id){
      rec_str=str.substr(str.find_first_of(" \t")+1);
      return;
    }
  }
  cout << "ERROR: no string with rec_id = " << rec_id << endl; exit(1);
}

// Allocates memmory for play_rec and last_motor_torques
void modelplayer::make_play_rec(){
  config_dim = model->get_config_dim();
  nmj = model->number_of_motor_joints();
  play_rec = new double [config_dim];
  last_motor_torques = new double [nmj];
}

// Prepares a (periodic cycle) trajectory of length n_t 
// (points per cycle) in periodic obj, generated by pgs.
// Computes relevant dynamic quantities, such as velocities
// and accelerations of model parts, but not of its
// generalized coordinates. Switches on torso penalties.
void modelplayer::prepare_per_traj_dyn(periodic& per, pergensetup* pgs, int n_t){
  per.record_trajectory(pgs,n_t);
  per.compute_dynrecs();
  per.compute_dynrec_ders();
  per.switch_torso_penalty(1,1);
}

// Mesures COT for given pgs over a cycle of size n_t.
// Optinally prints minimal contact force
// and maximal (non-slip) friction coefficient. 
double modelplayer::measure_cot(pergensetup* pgs, int n_t){
  periodic per (model);
  prepare_per_traj_dyn(per,pgs,n_t);
  double work = per.work_over_period();
  double weight = per.get_total_mass();
  double step_length = pgs->get_pergen()->get_step_length();
  double cot = work/(weight*step_length);
  //cout << "COT = " << cot << endl;
  if(contact_force_flag){
    double stat[2];
    per.get_contforce_stat(stat);
    cout << "min cfz = " << stat[0];
    cout << ", max mu = " << stat[1] << endl;
  }
  
  return cot;
}

// Writes pgs congig parameters to a string,
// in config string format.
void modelplayer::pergensu_config_string(pergensetup* pgs, string& str){
  str.clear();
  string fname = model->get_xmlfname();
  str += " xml_file " + fname;
  extvec pos, angles;
  double step_duration, TLh[3];
  pgs->get_config_params(pos,angles,step_duration,TLh);
  stringstream ss;
  for(int i=0;i<3;i++){ss << " " << pos.get_v(i);}
  str += " torso_pos" + ss.str();
  ss.str("");
  for(int i=0;i<3;i++){ss << " " << angles.get_v(i);}
  str += " torso_angles" + ss.str();
  ss.str("");
  ss << " step_duration " << step_duration << " period " << TLh[0] << " step_length " << TLh[1] << " step_height " << TLh[2];
  str += ss.str();
}

// Measures COT by sweeping through values of a parameter.
// Sweeps through n_val values from val0 to val1, for param_name
// parameter. Supported parameters: step_duration, period,
// step_length and step height.
void modelplayer::measure_cot_sweep(pergensetup* pgs, int n_t, string param_name, double val0, double val1, int n_val){

  pgssweeper sweeper (pgs, model);
  sweeper.sweep(param_name, val0, val1, n_val);
  while(sweeper.next()){
    pergensetup* pgs1 = sweeper.get_pgs();
    double cot = measure_cot(pgs1, n_t);
    double val = sweeper.get_val();
    cout << "val = " << val << " COT = " << cot << endl;
  }
}

// various simulation step updates, such as
// settin motor torques, followed by ODE
// simulation step.
void modelplayer::simulate_ode(){
  if(dynamics_from_simulation_flag){estimate_B();}
  //model->print(2);//exit(1);
  if(open_loop_flag){set_open_loop_torques();}
  if(position_control_flag){set_position_control_torques();}
  if(cpc_control_flag){set_cpc_torques();}
  if(record_traj_flag){add_traj_record();}
  if(torso_kicks_flag){kick_torso();}
  //torso_velocity(); // experim
  //print_array<double>(last_motor_torques,nmj); // temp
  //if(int(play_t/play_dt)%50==0){fall_check(.4);}
  get_vis()->simulate_odeworld(play_dt);
  play_t += play_dt;
}

// Sets drawing speedup factor f, which is the number 
// of simulation steps per one ds visualization.
void modelplayer::speedup_draw(int f){
  get_vis()->set_speedup(f);
}

// Sets initial configuration to be the state
// from pgs trajectory at time play_t, but with
// all velocities = 0.
void modelplayer::init_play_config(pergensetup* pgs){
  pgs->set_rec(play_rec,play_t);
  set_jangles_with_lik(play_rec);
  model->recompute_modelnodes();
  model->orient_odebodys();
}

void modelplayer::position_control_test(pergensetup* pgs, double t0){
  set_flag("position_control",true);
  setup_per_controller(pgs,t0);//exit(1);
  model->draw();
  unset_per_controller();
  set_flag("position_control",false);
}

// helper function to use in a controller setup.
// Constructs play_per, that makes trajectory and computes 
// forces over a period (cycle). Sets sim mode, sets initial
// state configuration (all velocities are zero).
void modelplayer::setup_per_controller(pergensetup* pgs, double t0){
  double T = pgs->get_period();
  int n_t = int(T/play_dt+.5);
  play_per = new periodic (model);
  prepare_per_traj_dyn(*play_per,pgs,n_t);
  play_per->compute_torques_over_period();
  play_per->compute_vel_traj();
  
  play_t = int(t0/play_dt+.5)*play_dt;
  //set_flag("manual_viewpoint",false);
  step_mode = 6;
  init_play_config(pgs);
}

void modelplayer::unset_per_controller(){
  delete play_per;
}

// Sets position control torques by: 
// reading targets state from play_per,
// reading current state from vis,
// computing linear feedback controls
// and setting motor torques.
void modelplayer::set_position_control_torques(){
  double k = 100;
  double k1 = -k, k2 = -2*sqrt(k);
  int an = 5;
  double** a = new_2d_array(an,nmj);
  double *q0 = a[0], *dq0 = a[1], *q = a[2], *dq = a[3];
  double *motor_torques = a[4];
  int tsi = int(play_t/play_dt+.5);
  play_per->get_motor_adas(tsi,q0,dq0);
  get_vis()->get_ode_motor_adas(q,dq);
  if(ghost_walking_flag){ghost->get_motor_adas(q,dq);}
  double *p = play_per->get_computed_torques(tsi);
  std::copy(p,p+nmj,motor_torques);
  double *x0[2] = {q0,dq0}, *x[2] = {q,dq};
  linear_feedback_control(motor_torques,x0,x,k1,k2);
  //arrayops ao (nmj); cout << ao.norm(p) << " " << ao.norm(motor_torques) << " " << ao.distance(motor_torques,p) << endl;
  //arrayops ao (nmj); cout << ao.norm(motor_torques)<<endl;
  set_ode_motor_torques(motor_torques);
  delete_2d_array(a,an);
}

// Computes linear feedback torques given target state x0,
// current state x, position and velocity gains k1 and k2.
void modelplayer::linear_feedback_control(double* torques, double** x0, double** x, double k1, double k2){
  double *q0 = x0[0], *dq0 = x0[1], *q = x[0], *dq = x[1];
  double* a1 = new double [2*nmj];
  double* a2 = a1 + nmj;
  arrayops ao (nmj);
  ao.assign(a1,q);
  ao.assign(a2,dq);
  //ao.assign_scalar(torques,0); // temp, test
  //ao.print(torques);ao.print(q0);ao.print(dq0);
  ao.modulus(ao.subtract(a1,q0),2*M_PI);
  ao.subtract(a2,dq0);
  //cout<<ao.l1_norm(a1)+ao.l1_norm(a2)*0<<endl;
  ao.add(ao.times(a1,k1),ao.times(a2,k2));
  ao.add(torques,a1);
  //ao.print(torques); cout<<endl;
  delete [] a1;
}

void modelplayer::save_last_motor_torques(const double* torques){
  arrayops ao (nmj);
  ao.assign(last_motor_torques,torques);
}

// Configuration path control (CPC) test.
void modelplayer::cpc_test(pergensetup* pgs, double t0){
  set_flag("dynamics_from_simulation",true);
  set_flag("cpc_control",true);
  setup_per_controller(pgs,t0);
  setup_cpc_controller();
  model->draw();
  unset_per_controller();
  set_flag("cpc_control",false);
  exit(1);
}

void modelplayer::setup_cpc_controller(){
  play_cpc = new cpccontroller (model);
  play_cpc->set_w_sg_nd_k0_kc_tauc(10,1,10,100*2,10,20); //(10*10,1,10,100*2,10,20)
  play_cpc->set_target_points_by_per(play_per);
  //play_cpc->load_tpset("traj.txt");
  play_cpc->setup_effdata();
  //play_cpc->set_flag("poscontrol",true); // experimental
  play_cpc->set_flag("goal_t0s",true); // experimental
}

// Sets CPC torques by:
// copying estimated B and current state to cpc controller,
// computing torques with CPC algorithm,
// setting motor torques. 
void modelplayer::set_cpc_torques(){
  //if(play_t>20){exit(1);}//cout<<play_t<<endl;
  //if(play_t<10){set_position_control_torques();return;}
  if(play_t<2*play_dt){set_position_control_torques();return;}

  play_cpc->set_B_by_player(B_transp);
  double** ders = tdertrack->get_ders();
  play_cpc->set_current_state(ders[0],ders[1]);
  if(ghost_walking_flag){set_ghost_cpc_state(ders);}
  double *motor_torques = new double [nmj];
  play_cpc->get_motor_torques(motor_torques);
  set_ode_motor_torques(motor_torques);
  //play_cpc->get_torques0(motor_torques); // experimental
  //arrayops ao (nmj); ao.times(motor_torques,0); // experimental
  delete [] motor_torques;
  //arrayops ao (nmj);cout<<ao.norm(last_motor_torques)<<endl;cout<<endl;//cin.ignore();
}

// Records trajectory produced by position control
// starting at t0 for duration of traj_duration.
void modelplayer::record_pos_control_traj(pergensetup* pgs, double t0, double traj_duration){
  traj_t_limit = t0 + traj_duration;
  set_flag("dynamics_from_simulation",true);
  set_flag("record_traj",true);
  position_control_test(pgs,t0);
  set_flag("record_traj",false);
}

// Adds a record to trajectory recording.
// When play_t >= traj_t_limit, trajectory
// gets saved and program exits.
void modelplayer::add_traj_record(){
  int rec_len = 2*config_dim+nmj;
  double* rec = new double [rec_len];
  double** ders = tdertrack->get_ders();
  arrayops ao (config_dim), ao1 (nmj);
  ao.assign(rec,ders[0]);
  ao.assign(rec+config_dim,ders[1]);
  ao1.assign(rec+2*config_dim,last_motor_torques);
  traj_recording.push_back(rec);
  if(play_t < traj_t_limit){return;}
  save_traj_recording("traj.txt",rec_len);
  exit(1);
}

// Saves recorded trajectory to a file.
void modelplayer::save_traj_recording(string fname, int rec_len){
  int size = traj_recording.size();
  double** traj = new double* [size];
  std::copy(traj_recording.begin(),traj_recording.end(),traj);
  int offset = 2;
  save_2d_array(traj+offset,size-offset,rec_len,"traj.txt",false);
  delete [] traj;
}

// Sets right m columns of B transposed to unity.
// Note, n > m;
void set_Bt_unity(double** Bt, int m, int n){
  arrayops ao (n);
  double** p = Bt;
  for(int i=0;i<m;i++){
    ao.assign_scalar(*p,0);
    (*p++)[n-m+i] = 1;
  }
}

// Estimates B from T and U, (see CPC paper notations).
void B_from_T_U(double** Bt, MatrixXd& T, MatrixXd& U){
  int nmj = T.rows(), config_dim = U.rows();
  if(T.cols()==0){set_Bt_unity(Bt,nmj,config_dim);return;}
  T.conservativeResize(nmj+1,T.cols());
  T.row(nmj).setOnes();
  ColPivHouseholderQR<MatrixXd> decomp(T.transpose());
  MatrixXd B (config_dim,nmj+1);
  for(int i=0;i<config_dim;i++){
    B.row(i) = decomp.solve(U.row(i).transpose());
  }
  //cout<<B<<endl;
  for(int i=0;i<nmj;i++){
    VectorXd::Map(Bt[i],config_dim) = B.col(i);
  }
}

// Runs m simulations from the current state,
// estimted B by regression.
void modelplayer::estimate_B_by_regression(){
  int m = config_dim + 10*2;//*2; // experim
  // m=0; // for unity Bt
  double* config = new double [config_dim];
  double* torques = new double [nmj];
  arrayops ao (nmj);
  visualizer* vis = get_vis();
  vis->get_ode_config(config);
  tdertrack->push_config(config);
  odestate state (model);
  state.save();
  MatrixXd T(nmj,m), U(config_dim,m);
  double del_torque = .2*2;//*2; // experim
  for(int i=0;i<m;i++){
    configtimedertrack tdertrack1 (tdertrack);
    ao.assign(torques,last_motor_torques);
    for(int j=0;j<nmj;j++){
      torques[j] += float(rand()%100-50)/50.*del_torque;
    }
    vis->set_ode_motor_torques(torques);
    vis->simulate_odeworld(play_dt);
    vis->get_ode_config(config);
    tdertrack1.push_config(config);
    double** ders = tdertrack1.get_ders();
    T.col(i) = Map<VectorXd> (torques,nmj);
    U.col(i) = Map<VectorXd> (ders[2],config_dim);
    state.load();
  }
  //int a[]={0,1,5};for(int i=0;i<3;i++){U.row(a[i])*=0;} // experim, probably to remove, it is not needed
  delete [] config;
  delete [] torques;
  B_from_T_U(B_transp,T,U);
}

// Torso perturbation by impuls forces ("kicks").
void modelplayer::kick_torso(){
  dBodyID odebody = get_torso_odebody();
  //dBodyID odebody = (*model->get_odeparts())[12]->get_odebody();
  int k = 500;
  if(int(play_t/play_dt)%k!=(k-1)){return;}
  float th = 2*3.1416*float(rand()%100)/100.;
  float kick = float(rand()%100)/100.;
  kick=4.;th=M_PI*(rand()%2);
  th=0;
  //cout<<kick<<" "<<f<<endl;
  float dvx = kick*cos(th), dvz = kick*sin(th);
  cout<<dvx<<endl;
  //const dReal* v = dBodyGetLinearVel(body);
  //dBodySetLinearVel(body,v[0]+dvx,v[1],v[2]+dvz);

  double f[] = {dvx/play_dt,0,dvz/play_dt};
  get_vis()->add_force(odebody,f);
  //dBodyAddForce(odebody,f[0],f[1],f[2]);
  //dBodyAddForce(odebody,dvx/play_dt,0,dvz/play_dt);
  //dBodyAddTorque(odebody,0,dvx/play_dt,0);
}

void modelplayer::set_default_flags(){
  manual_viewpoint_flag = true, contact_force_flag = false, open_loop_flag = false, position_control_flag = false, dynamics_from_simulation_flag = false, cpc_control_flag = false, record_traj_flag = false, torso_kicks_flag = false, ghost_walking_flag = false, clip_torque_flag = false;
}

// Loads and draws model as defined in xml file.
void modelplayer::load_draw_model(string fname){
  load_model(fname);
  test(7);
}

// Records trajectory as generated by periodic.
// Trajectory need not be physically realizable.
void modelplayer::record_per_traj(pergensetup* pgs){
  double T = pgs->get_period();
  int n_t = int(T/play_dt+.5);
  int rec_len = 2*config_dim+nmj;
  double** traj = new_2d_array(n_t,rec_len);
  periodic* per = new periodic (model);
  prepare_per_traj_dyn(*per,pgs,n_t);
  per->get_complete_traj(traj);
  save_2d_array(traj,n_t,rec_len,"traj.txt",false);
  delete per;
  delete_2d_array(traj,n_t);
}

// Records n_val trajectories by sweeping over
// a parameter param_name, from val0 to val1.
void modelplayer::record_per_traj_sweep(pergensetup* pgs, string param_name, double val0, double val1, int n_val){
  double T = pgs->get_period();
  int n_t = int(T/play_dt+.5);
  int rec_len = 2*config_dim+nmj;
  double** traj = new_2d_array(n_t,rec_len);

  pgssweeper sweeper (pgs, model);
  sweeper.sweep(param_name, val0, val1, n_val);
  bool flag = false;
  while(sweeper.next()){
    pergensetup* pgs1 = sweeper.get_pgs();
    //cout<<"val = "<<sweeper.get_val()<<endl;
    periodic* per = new periodic (model);
    prepare_per_traj_dyn(*per,pgs1,n_t);
    per->get_complete_traj(traj);
    save_2d_array(traj,n_t,rec_len,"traj.txt",flag);
    delete per;
    if(!flag){flag = true;}
  }

  delete_2d_array(traj,n_t);
}

// Outputs torso velocity.
void modelplayer::torso_velocity(){
  dBodyID odebody = get_torso_odebody();
  const dReal* vel = dBodyGetLinearVel(odebody);
  print_array<const dReal>(vel,3,"torso vel: ");
}

dBodyID modelplayer::get_torso_odebody(){
  return get_vis()->get_torso_opart()->get_odebody();
}

// Detects model fall: if torso height gets < hc, program aborts.
void modelplayer::fall_check(double hc){
  const dReal* pos = dBodyGetPosition(get_torso_odebody());
  print_array<const dReal>(pos,2);
  if(pos[2] < hc){cout << "Fall at t = " << play_t << endl; exit(1);}
}

// Test for navigating various kinds of uneven ground. 
void modelplayer::uneven_ground_test(){
  float f = 1;// .2;
  double l = 1*f;
  int n = 10/f;
  hfield = new heightfield (n,5*1/f+2,l);
  //hfield->random_field(0.01,.3+.4,false);
  //hfield->random_field(0.01,.3,false);
  //hfield->random_field(0.01,.01,false);
  //hfield->slope_field(0,6.5);
  //hfield->tan_field(1.5,4);
  //hfield->tanh_field(5);
  //hfield->gauss_field(4);
  hfield->ridge_field(0,2.0);
  /* hfield->slope_field(0,4);
     hfield->ripple_field(0.01,.4*f,true); */
  //hfield->random_field(0.01,.2,true);
  //hfield->make_geom(n*l/2+3,0,get_vis());
  visualizer* vis = get_vis();
  vis->push_geom(hfield->make_geom(n*l/2+3-1.5,0,vis->get_trimeshman()));

  //return;
  set_flag("ghost_walking", true);
  ghost = new ghostmodel (vis, hfield, play_dt);
  
  extvec eas (0,0.0,0);
  ghost->set_surf_rot(eas);

  set_torque_limit(10);
}

// Sets ghost state as the current state in CPC controller.
void modelplayer::set_ghost_cpc_state(double** ders){
  double** a = new_2d_array(2,config_dim);
  double *q = a[0], *dq = a[1];
  std::copy(ders[0],ders[0]+6,q);
  std::copy(ders[1],ders[1]+6,dq);
  ghost->get_motor_adas(q+6,dq+6);
  play_cpc->set_current_state(q,dq);
  delete_2d_array(a,2);
}

// Outputs maximal motor torque component.
void modelplayer::max_motor_torque(){
  float max_torque = 0;
  double* p = last_motor_torques;
  for(int i=0;i<nmj;i++){
    double torque = fabs(*p++);
    if(torque > max_torque){max_torque = torque;}
  }
  cout << "max torque = " << max_torque << endl;
}

void modelplayer::set_ode_motor_torques(double* motor_torques){
  if(clip_torque_flag){clip_torque(motor_torques);}
  get_vis()->set_ode_motor_torques(motor_torques);
  save_last_motor_torques(motor_torques);
  //max_motor_torque(); // experim
}

void modelplayer::set_torque_limit(double torque){
  torque_limit = torque;
  set_flag("clip_torque",true);
}

// Clips torque components exceeding torque_limit.
void modelplayer::clip_torque(double* motor_torques){
  double* p = motor_torques;
  for(int i=0;i<nmj;i++){
    double torque = *p;
    if(fabs(torque) > torque_limit){
      *p = torque_limit*torque/fabs(torque);
    }
    p++;
  }
}

void modelplayer::test_lik_solvers(){
  model->get_lik()->solver_test(10000);
}

// Shifts camera view point.
void modelplayer::shift_view(double x, double y, double z){
  get_vis()->get_view()->shift_cam(x,y,z);
}

// Checks if model is loaded, aborts othersize.
void modelplayer::check_model_loaded(){
  if(!model->if_loaded()){cout<<"ERROR: model not loaded"<<endl;exit(1);}
}


/*

//euler angles check, from player:
affine* A = get_vis()->get_torso_opart()->get_mnode()->get_A_ground(); double as[3]; euler_angles_from_affine(*A, as); print_array<double>(as,3); exit(1);

*/

