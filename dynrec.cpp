#include "dynrec.h"


void dynpart::setup(map<modelnode*,int>& mnode_id_map){
  modelnode* pmnode;
  mnode = opart->get_mnode();
  pmnode = mnode->get_parent();
  id = mnode_id_map[mnode];
  parent_id = (pmnode)? mnode_id_map[pmnode] : -1;
  recompute();
  set_inertial_params();
}

void dynpart::print(){
  cout << "--- dynpart ---" << endl;
  cout << "name: " << opart->get_part_name() << endl;
  cout << "id = " << id << endl;
  cout << "parent id = " << parent_id << endl;
  cout << "joint_pos: ";
  joint_pos.print();
  cout << "com_pos: ";
  com_pos.print();
  if(foot_flag){
    cout << "foot_pos: ";
    foot_pos.print();
  }
}

void dynpart::set_joint_pos(){
  modeljoint* joint = mnode->get_joint();
  if(joint){
    joint->get_A_ground()->get_translation(joint_pos);
  } else {
    mnode->get_A_ground()->get_translation(joint_pos);
  }
  //joint_pos.print();
}

void dynpart::set_com_pos(){
  opart->get_com_pos(com_pos);
}

void dynpart::recompute(){
  set_joint_pos();
  set_com_pos();
  if(foot_flag){set_foot_pos();}
}

const affine* dynpart::get_A_ground(){
  return mnode->get_A_ground();
}

affine* dynpart::get_inertia_tensor(){
  return &inertia;
}

void dynpart::set_inertial_params(){
  dMass m;
  dBodyGetMass(opart->get_odebody(),&m);
  mass = m.mass;
  inertia.set_rotation(m.I);
  //cout<<mass<<endl;inertia.print();
}

void dynpart::setup_foot(const set<modelnode*>& foot_set){
  if(foot_set.count(mnode)){
    foot_flag = true;
    set_foot_pos();
  } else {foot_flag = false;}
}

void dynpart::set_foot_pos(){
  opart->get_foot_pos(foot_pos);
}

void dynpart::get_joint_zaxis(extvec& axis){
  modeljoint* joint = mnode->get_joint();
  if(joint){
    affine* A = joint->get_A_ground();
    double *p = A->get_data() + 8;
    axis.set(p);
  } else {
    axis.set_zeros();
  }
}


dynrecord::dynrecord(int n_, int nf_){
  n = n_;
  nf = nf_;
  pos = new extvec [n];
  jpos = new extvec [n];
  vel = new extvec [n];
  mom = new extvec [n];
  mom_rate = new extvec [n];
  acc = new extvec [n];
  ust = new extvec [n];
  ang_vel = new extvec [n];
  ang_mom = new extvec [n];
  ang_mom_rate = new extvec [n];
  fpos = new extvec [nf];
  jzaxis = new extvec [n];
  rot = new affine [n];
  contacts = new bool [nf];
}

dynrecord::~dynrecord(){
  delete [] pos;
  delete [] jpos;
  delete [] vel;
  delete [] mom;
  delete [] mom_rate;
  delete [] acc;
  delete [] ust;
  delete [] ang_vel;
  delete [] ang_mom;
  delete [] ang_mom_rate;
  delete [] fpos;
  delete [] jzaxis;
  delete [] rot;
  delete [] contacts;
}

void dynrecord::initialize(vector<dynpart*>& dynparts, double rcap){
  vector<dynpart*>::iterator it = dynparts.begin();
  int fi = 0;
  for(int i=0;i<n;i++){
    dynpart* dpart = (*it);
    //pos[i].copy(*dpart->get_com_pos());
    //jpos[i].copy(*dpart->get_joint_pos());
    pos[i] = (*dpart->get_com_pos());
    jpos[i] = (*dpart->get_joint_pos());
    const affine* A = dpart->get_A_ground();
    double x = (A->get_a(2,1)-A->get_a(1,2))/2;
    double y = (A->get_a(0,2)-A->get_a(2,0))/2;
    double z = (A->get_a(1,0)-A->get_a(0,1))/2;
    ust[i].set(x,y,z);
    rot[i].set_rotation(*A);
    if(dpart->if_foot()){
      //fpos[fi].copy(*dpart->get_foot_pos());
      fpos[fi] = (*dpart->get_foot_pos());
      contacts[fi] = (fpos[fi].get_data()[2] < rcap+1e-4);
      fi++;
    }
    dpart->get_joint_zaxis(jzaxis[i]);
    it++;
  }
}

void dynrecord::print(int id){
  cout << "--- dynrec ---" << endl;
  cout << "pos: "; pos[id].print();
  cout << "jpos: "; jpos[id].print();
  cout << "vel: "; vel[id].print();
  cout << "mom: "; mom[id].print();
  cout << "mom_rate: "; mom_rate[id].print();
  cout << "acc: "; acc[id].print();
  cout << "ust: "; ust[id].print();
  cout << "ang_vel: "; ang_vel[id].print();
  cout << "ang_mom: "; ang_mom[id].print();
  cout << "ang_mom_rate: "; ang_mom_rate[id].print();
  cout << "rot: " << endl; rot[id].print();
}

void dynrecord::compute_ders(int stage, const dynrecord* prev_rec, const dynrecord* next_rec, double dt, periodic* per){
  switch (stage) {
  case 0: {
    compute_ders(vel,prev_rec->pos,next_rec->pos,dt);
    compute_ders(ang_vel,prev_rec->ust,next_rec->ust,dt);
    compute_mom(per->get_masses());
    compute_ang_mom(per);
  } break;    
  case 1: {
    compute_ders(mom_rate,prev_rec->mom,next_rec->mom,dt);
    compute_ders(acc,prev_rec->vel,next_rec->vel,dt); // probably acc wont be needed
    compute_ders(ang_mom_rate,prev_rec->ang_mom,next_rec->ang_mom,dt);
  } break;    
  }
}

void dynrecord::compute_ders(extvec *p_der, const extvec *p_func_prev, const extvec *p_func_next, double dt){
  for(int i=0;i<n;i++){
    //p_der->copy(*p_func_next);
    *p_der = (*p_func_next);
    p_der->subtract(*p_func_prev);
    p_der->times(1./(2*dt));
    
    p_der++;
    p_func_prev++;
    p_func_next++;
  }
}

void dynrecord::compute_ang_mom(periodic* per){
  for(int i=0;i<n;i++){
    affine* A_I = per->get_dynpart(i)->get_inertia_tensor();
    extvec v, u;
    affine rot_tr;
    rot_tr.copy_transposed(rot[i]);
    rot_tr.mult(ang_vel[i],v);
    A_I->mult(v,u);
    rot[i].mult(u,ang_mom[i]);
  }
}

void dynrecord::compute_mom(const double* masses){
  for(int i=0;i<n;i++){
    double mass = masses[i];
    //mom[i].copy(vel[i]);
    mom[i] = vel[i];
    mom[i].times(mass);
  }
}

void dynrecord::set_forcetorque_system(SpMat& B, VectorXd& f, const int* parentis, const double* masses){
  for(int i=0;i<n;i++){
    int pi = parentis[i];
    ftsys_forces(i,pi,B,f);
    ftsys_torques(i,pi,B,f);
    ftsys_gravity(i,f,masses);
  }
  B.makeCompressed();
}

void ftsys_unit_elems(int i, int pi, SpMat& B){
  for(int j=0;j<3;j++){
    int k = 3*i+j, k1 = 3*pi+j;
    B.insert(k,k) = 1;
    if(pi >= 0){B.insert(k1,k) = -1;}
  }
}

// effect of force applied at joint j on part i
void ftsys_cross_elems(int i, int j, extvec& r, SpMat& B, int n){
  int k = 3*(n+i), k1 = 3*j;
  double *p = r.get_data();
  int dk[3];
  for(int l=0;l<3;l++){
    for(int l1=0;l1<3;l1++){dk[l1] = (l+l1)%3;}
    B.insert(k+dk[0],k1+dk[1]) = -p[dk[2]];   
    B.insert(k+dk[1],k1+dk[0]) = p[dk[2]];   
  }
}

void ftsys_cross_elems(int i, int pi, extvec& ipos, extvec& jpos, extvec* ppos, SpMat& B, int n){
  extvec r (jpos);
  r.subtract(ipos);
  ftsys_cross_elems(i,i,r,B,n);
  //r.copy(*ppos);
  r = (*ppos);
  r.subtract(jpos);
  ftsys_cross_elems(pi,i,r,B,n);
}

void dynrecord::ftsys_forces(int i, int pi, SpMat& B, VectorXd& f){
  ftsys_unit_elems(i,pi,B);
  int k0 = 3*i;
  mom_rate[i].get_components(f(k0),f(k0+1),f(k0+2));
}

void dynrecord::ftsys_torques(int i, int pi, SpMat& B, VectorXd& f){
  int pi1 = (pi < 0)? pi : pi+n;
  ftsys_unit_elems(i+n,pi1,B);
  if (pi >= 0) {ftsys_cross_elems(i,pi,pos[i],jpos[i],&(pos[pi]),B,n);}
  int k0 = 3*(n+i);
  ang_mom_rate[i].get_components(f(k0),f(k0+1),f(k0+2));  
}

void dynrecord::ftsys_gravity(int i, VectorXd& f, const double* ms){
  const double g = 1;
  //const double g = 10;
  f(3*i+2) += ms[i]*g;
}

int dynrecord::get_ncontacts(){
  int s = 0;
  bool *p = contacts;
  for(int i=0;i<nf;i++){s += *p++;}
  return s;
}

void dynrecord::set_forcetorque_system_contacts(SpMat& B, const int* footis){
  set_forcetorque_system_contacts(B,footis,true);
}

void dynrecord::set_forcetorque_system_contacts(SpMat& B, const int* footis, bool contact_feet_flag){
  int ci = 0;
  for(int fi=0;fi<nf;fi++){
    if(contact_feet_flag){
      if(!contacts[fi]){continue;}
    }
    int i = footis[fi];
    ftsys_contact_forces(i,ci,B);
    ftsys_contact_torques(fi,i,ci,B);
    ci++;
  }
  B.makeCompressed();
}

// for i and ci see ftsys_contact_torques comment
void dynrecord::ftsys_contact_forces(int i, int ci, SpMat& B){
  for(int j=0;j<3;j++){
    int  k = 3*i+j, k1 = 3*(2*n+ci)+j;
    B.insert(k,k1) = 1;
  }
}

// i - foot part id
// fi - foot id
// ci - contact id
void dynrecord::ftsys_contact_torques(int fi, int i, int ci, SpMat& B){
  extvec r (fpos[fi]);
  r.subtract(pos[i]);
  ftsys_cross_elems(i,2*n+ci,r,B,n);
}

