#include "odestate.h"


void odebodydata::read(){
  const dReal* odebody_pos = dBodyGetPosition(odebody);
  const dReal* odebody_vel = dBodyGetLinearVel(odebody);
  const dReal* odebody_ang_vel = dBodyGetAngularVel(odebody);
  const dReal* odebody_quat = dBodyGetQuaternion(odebody);
  dReal *p = pos, *p1 = vel, *p2 = ang_vel, *p3 = quat;
  for(int i=0;i<4;i++){
    *p++ = *odebody_pos++;
    *p1++ = *odebody_vel++;
    *p2++ = *odebody_ang_vel++;
    *p3++ = *odebody_quat++;
  }
}

void odebodydata::write(){
  //float e = float(rand()%100-50)/100.;
  dBodySetPosition(odebody,pos[0],pos[1],pos[2]);
  dBodySetLinearVel(odebody,vel[0],vel[1],vel[2]);
  dBodySetAngularVel(odebody,ang_vel[0],ang_vel[1],ang_vel[2]);
  const dQuaternion odebody_quat = {quat[0],quat[1],quat[2],quat[3]};
  dBodySetQuaternion(odebody,odebody_quat);
}

void odebodydata::print(){
  print_array<dReal>(pos,4,"pos: ");
  print_array<dReal>(vel,4,"vel: ");
  print_array<dReal>(ang_vel,4,"ang_vel: ");
  print_array<dReal>(quat,4,"quat: ");
}


odestate::odestate(const kinematicmodel* model){
  const vector<odepart*>* odeparts = model->get_odeparts();
  vector<odepart*>::const_iterator it = odeparts->begin();
  for(;it!=odeparts->end();it++){
    odebodydata* obody = new odebodydata ((*it)->get_odebody());
    obodys.push_back(obody);
  }
}

odestate::~odestate(){
  list<odebodydata*>::iterator it = obodys.begin();
  for(;it!=obodys.end();it++){delete *it;}
}

void odestate::save(){
  list<odebodydata*>::iterator it = obodys.begin();
  for(;it!=obodys.end();it++){(*it)->read();}
}

void odestate::load(){
  list<odebodydata*>::iterator it = obodys.begin();
  for(;it!=obodys.end();it++){(*it)->write();}
}

void odestate::print(){
  cout << "----- ode state -----" << endl;
  list<odebodydata*>::iterator it = obodys.begin();
  int i = 0;
  for(;it!=obodys.end();it++){
    cout << "obody " << i++ << ":" << endl;
    (*it)->print();
  }
}


configtimedertrack::configtimedertrack(int max_der_, int config_dim_, double dt_){
  construct(max_der_,config_dim_,dt_);
}

configtimedertrack::configtimedertrack(const configtimedertrack* tdertrack){
  const configtimedertrack* tdt = tdertrack;
  construct(tdt->max_der,tdt->config_dim,tdt->dt);
  copy(tdt);
}

void configtimedertrack::construct(int max_der_, int config_dim_, double dt_){
  max_der = max_der_;
  config_dim = config_dim_;
  dt = dt_;
  ders = new_2d_array(max_der+1,config_dim);
  new_der = new double [config_dim];
  old_der = new double [config_dim];
  ao.set_n(config_dim);
  for(int i=0;i<=max_der;i++){ao.assign_scalar(ders[i],0);}
}

configtimedertrack::~configtimedertrack(){
  delete_2d_array(ders,max_der+1);
  delete [] new_der;
  delete [] old_der;
}

void configtimedertrack::push_config(double* config){
  ao.assign(old_der,ders[0]);
  ao.assign(ders[0],config);
  ao.assign(new_der,config);
  for(int i=0;i<max_der;i++){compute_dern(i+1);}
}

// use it for n > 0
void configtimedertrack::compute_dern(int n){
  ao.subtract(new_der,old_der);
  ao.times(new_der,1./dt);
  ao.assign(old_der,ders[n]);
  ao.assign(ders[n],new_der);
}

void configtimedertrack::print(){
  cout << "----- configuration time ders -----" << endl;
  for(int i=0;i<=max_der;i++){
    cout << "der " << i << ": ";
    print_array<double>(ders[i],config_dim);
  }
}

void configtimedertrack::copy(const configtimedertrack* tdertrack){
  const configtimedertrack* tdt = tdertrack;
  for(int i=0;i<max_der;i++){ao.assign(ders[i],tdt->ders[i]);}
  ao.assign(new_der,tdt->new_der);
  ao.assign(old_der,tdt->old_der);
}

void configtimedertrack::del_second_der(configtimedertrack* tdertrack, double* del){
  ao.assign(del,ders[2]);
  ao.subtract(del,tdertrack->ders[2]);
}

