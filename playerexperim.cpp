#include "lik.h"
#include "pergen.h"
#include "player.h"
#include "periodic.h"
#include "odestate.h"
#include "ghost.h"

// quadruped spinning with foot tips fixed
void modelplayer::test0(){
  vector<double*>* joint_values = model->get_joint_values();
  for(int i=0;i<(int)joint_values->size();i++){
    if(i==5){*(*joint_values)[i]+=.02;}
  }
  const liksolver* lik = model->get_lik();
  for(int i=0;i<4;i++){lik->place_limb(i,.25*cos(i*M_PI/2),.25*sin(i*M_PI/2),1);}
}

// quadruped twisting
void modelplayer::test1(){
  vector<double*>* joint_values = model->get_joint_values();
  double delz = -.1;
  play_t += .05;
  *(*joint_values)[2] = delz;
  *(*joint_values)[5] = .5*sin(play_t);
model->recompute_modelnodes();
  const liksolver* lik = model->get_lik();
  for(int i=0;i<4;i++){
    double a = (i+.5)*M_PI/2;
    lik->place_limb(i,.4*cos(a),.4*sin(a),.2+delz);
  }
}

// hexapod
void modelplayer::test2(){
  if(play_t==0){
    play_pgs = new pergensetup (6);
    extvec torso_pos (0,0,-0.1);
    extvec euler_angles (0,0,0*M_PI/2);
    extvec orientation[] = {torso_pos,euler_angles};
    setup_pergen(*play_pgs,orientation,.5);
    play_pgs->set_TLh(.5,.8,.1);
  }
  play_pgs->set_rec(play_rec,play_t);
  set_jangles_with_lik(play_rec);
  play_t += play_dt;
}

// quadruped 2-leg walking
void modelplayer::test3(){
  if(play_t==0){
    play_pgs = new pergensetup (4);
    //extvec torso_pos (0,0,-.5);
    //extvec euler_angles (0,0,0);
    extvec torso_pos (0,0,.35);
    extvec euler_angles (0,M_PI/2,0*M_PI/2);
    extvec orientation[] = {torso_pos,euler_angles};
    setup_pergen(*play_pgs,orientation,1);
    extvec delpos0 (0,0,.9);
    for(int i=1;i<4;i+=2){
      play_pgs->get_pergen()->change_pos0(i,delpos0);
    }
    play_pgs->set_TLh(.5,.6,.1);
  }
  play_pgs->set_rec(play_rec,play_t);
  set_jangles_with_lik(play_rec);
  play_t += play_dt;
}

// quadruped
void modelplayer::test4(){
  if(play_t==0){
    play_pgs = new pergensetup (4);
    //extvec torso_pos (0,0,-.1);
    //extvec euler_angles (0,0,0);
    extvec torso_pos (0,0,-.07);
    extvec euler_angles (0,0,-.5*M_PI/2);
    extvec orientation[] = {torso_pos,euler_angles};
    setup_pergen(*play_pgs,orientation,.5);
    play_pgs->set_TLh(.5,.7,.1);
  }
  play_pgs->set_rec(play_rec,play_t);
  set_jangles_with_lik(play_rec);
  play_t += play_dt;
}

// open loop controller
void modelplayer::open_loop_test(pergensetup* pgs, double t0){
  set_flag("open_loop",true);
  setup_per_controller(pgs,t0);
  model->draw();
  unset_per_controller();
  set_flag("open_loop",false);
}

// blindly applies computed torques, no feedback
void modelplayer::set_open_loop_torques(){
  //cout<<play_t<<endl;
  int i = int(play_t/play_dt+.5);
  double* motor_torques = play_per->get_computed_torques(i);
  set_ode_motor_torques(motor_torques);
}

void modelplayer::estimate_B(){
  estimate_B_by_regression();return; // temporary/experimental
  //if(play_t>30){exit(1);}//cout<<play_t<<endl;
  double* config = new double [config_dim];
  double* torques = new double [nmj];
  arrayops ao (nmj), ao1 (config_dim);
  visualizer* vis = get_vis();
  vis->get_ode_config(config);
  tdertrack->push_config(config);
  odestate state (model);
  state.save();
  double del_torque = 1;
  //MatrixXd B (config_dim,nmj);
  configtimedertrack tdertrack0 (tdertrack);
  for(int i=0;i<=nmj;i++){
    configtimedertrack tdertrack1 (tdertrack);
    ao.assign(torques,last_motor_torques);
    if(i){torques[i-1] += del_torque;}
    vis->set_ode_motor_torques(torques);
    vis->simulate_odeworld(play_dt);
    vis->get_ode_config(config);
    if(i){
      tdertrack1.push_config(config);
      //tdertrack1.del_second_der(&tdertrack0,config);
      tdertrack1.del_second_der(&tdertrack0,B_transp[i-1]);
      ao1.times(B_transp[i-1],1./del_torque);
      //B.col(i-1) = Map<VectorXd> (B_transp[i-1],config_dim);
    } else {tdertrack0.push_config(config);}
    state.load();
  }
  state.load();
  //tdertrack->print();if(rand()%5==0){exit(1);}
  //print_array<double>(config,config_dim);exit(1);
  delete [] config;
  delete [] torques;
  //state.print();
  //B /= del_torque;
  //if(play_t>.05){cout<<B.leftCols(6)<<endl;exit(1);}
}
