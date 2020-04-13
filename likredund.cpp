#include "lik.h"
#include "likredund.h"

// Solves 4-link limb with hinge axis y-x-x-x, (weaver configuration).
// Foot position pos_limb is in the hip joint's frame.
// To change knee bend (inward/outward) flip sign of s1.
bool limb_solver_yxxx(const double* constrs, double* jangles, const double* ls, int ysign, bool bend){
  extvec pos_limb, joint_angles;
  pos_limb.set(constrs), joint_angles.set(jangles);
  //double l0=ls[0], l1=ls[1], l2=ls[2];
  double l0=ls[0], l1=ls[1], l2= ls[2]+ls[3];
  //double l00=ls[0], l0=ls[1], l1=ls[2], l2= ls[3]+ls[4];

  int s0 = ysign;
  int s1 = 2*int(bend)-1;

  extvec pos0 (0,0,s0*l0);
  //extvec pos0 (0,l00,s0*l0);
  extvec pos1 (pos_limb);
  pos1.subtract(pos0);
  double l = pos1.norm();
  if(l1+l2-l<0){
    if(ignore_reach_flag){l=l1+l2;}
    else {cout<<"LIK ERROR: limb position is unreachable"<<endl;return false;}
  }

  double x1, y1, z1; // rotated frame, x->x1, y->z1, z->-y1 (?)  
  pos_limb.get_components(x1,y1,z1);
  double c = (z1-s0*l0)/l;
  double phi = atan2(x1,y1);
  double theta = acos(c) + (1-s0)*M_PI/2;
  //if (fabs(phi)>M_PI/2) {phi += M_PI; theta *= -1;}

  mod_twopi(phi);
  mod_twopi(theta);

  double ll = l*l;
  double del = l2*l2-l1*l1;
  double beta = s1*s0*acos((ll-del)/(2*l1*l));
  double gamma = s1*s0*acos((ll+del)/(2*l2*l));
  joint_angles.set(-phi,-theta+beta,-(beta+gamma));
  //cout<<"ll="<<ll<<" del="<<del<<" l1="<<l1<<" l="<<l<<endl;
  //cout<<"beta="<<beta<<" theta="<<theta<<" gamma="<<gamma<<endl;

  //joint_angles.get_components(jangles+1);
  //jangles[0] = 0;
  //jangles[4] = 0;
  joint_angles.get_components(jangles);
  jangles[3] = 0;
  //for(int i=0;i<5;i++){jangles[i]=0;}
  //print_array(jangles,5);if(limbi==5){cin.ignore();}
  //joint_angles.print();
  return true;
}

// Solves 5-link limb with hinge axis z-y-x-x-x, (weaver configuration).
// Foot position pos_limb is in the hip joint's frame.
// To change knee bend (inward/outward) flip sign of s1.
bool limb_solver_zyxxx(const double* constrs, double* jangles, const double* ls, int ysign, bool bend){
  extvec pos_limb, joint_angles;
  pos_limb.set(constrs), joint_angles.set(jangles);
  //double l0=ls[0], l1=ls[1], l2=ls[2];
  //double l0=ls[0], l1=ls[1], l2= ls[2]+ls[3];
  double l00=ls[0], l0=ls[1], l1=ls[2], l2= ls[3]+ls[4];

  int s0 = ysign;
  int s1 = 2*int(bend)-1;
  //s1*=-1;

  //extvec pos0 (0,0,s0*l0);
  extvec pos0 (0,s0*l0,-l00);
  extvec pos1 (pos_limb); //pos_limb.print();exit(1);
  pos1.subtract(pos0);
  double l = pos1.norm();
  if(l1+l2-l<0){
    if(ignore_reach_flag){l=l1+l2;}
    else {cout<<"LIK ERROR: limb position is unreachable"<<endl;return false;}
  }

  double x1, y1, z1; // rotated frame, x->x1, y->z1, z->-y1 (?)  
  /*pos_limb.get_components(x1,y1,z1);
  double c = (z1-s0*l0)/l;
  double phi = atan2(x1,y1);*/
  pos_limb.get_components(x1,y1,z1);
  //z1 *= -1;z1 -= l00;
  z1 = -(z1+l00);
  double c = (y1-s0*l0)/l;
  double phi = atan2(x1,z1);
  double theta = acos(c) + (1-s0)*M_PI/2;
  //double theta = acos(c) + s0*M_PI/2;
  //if (fabs(phi)>M_PI/2) {phi += M_PI; theta *= -1;}

  mod_twopi(phi);
  mod_twopi(theta);

  double ll = l*l;
  double del = l2*l2-l1*l1;
  double beta = s1*s0*acos((ll-del)/(2*l1*l));
  double gamma = s1*s0*acos((ll+del)/(2*l2*l));
  joint_angles.set(-phi,-theta+beta,-(beta+gamma));
  //joint_angles.set(-phi,theta-beta,(beta+gamma));
  //cout<<"ll="<<ll<<" del="<<del<<" l1="<<l1<<" l="<<l<<endl;
  //cout<<"beta="<<beta<<" theta="<<theta<<" gamma="<<gamma<<endl;

  joint_angles.get_components(jangles+1);
  jangles[0] = 0;
  jangles[4] = 0;
  //joint_angles.get_components(jangles);jangles[3] = 0;
  //for(int i=0;i<5;i++){jangles[i]=0;}
  //print_array(jangles,5);if(limbi==5){cin.ignore();}
  //joint_angles.print();
  return true;
}

// Solves 5-link limb with hinge axis z-y-x-x-x, (weaver configuration).
// Foot position pos_limb is in the hip joint's frame.
// To change knee bend (inward/outward) flip sign of s1.
bool limb_solver_zyxxx1(const double* constrs, double* jangles, const double* ls, int ysign, bool bend){
  extvec pos_limb, joint_angles;
  pos_limb.set(constrs), joint_angles.set(jangles);
  //double l0=ls[0], l1=ls[1], l2=ls[2];
  //double l0=ls[0], l1=ls[1], l2= ls[2]+ls[3];
  //double l00=ls[0], l0=ls[1], l1=ls[2], l2= ls[3]+ls[4];
  double l0 = ls[0], l4 = ls[4];

  int s0 = ysign;

  //print_array(constrs,5);exit(1);
  double q[5];
  // step 1
  extvec n3, r0 (0,0,-l0), r4 (pos_limb);
  direction_to_vec(constrs+3,n3);
  n3.times(-1);
  //n.print();exit(1);
  // step 2
  extvec v (r4);
  v.subtract(r0);
  v.cross(n3);
  v.cross(r0);
  v.times(-1);
  // step 3
  double dir_n1[2];
  vec_to_direction(v,dir_n1);
  //q[0] = dir_n1[1]-M_PI/2;
  q[0] = dir_n1[1]-s0*M_PI/2;
  //print_array(dir_n1,2);
  affine S, A;
  set_xyz_rotation(S,2,q[0]);
  //set_xyz_rotation(A,0,-s0*M_PI/2);
  set_xyz_rotation(A,0,-1*M_PI/2);
  S.translate(r0);
  S.mult(A);
  //S.print();//v.print();exit(1);
  // step 4
  extvec r3 (r4), r3_p; // r3_p = "r3 prime"
  v.copy(n3);
  v.times(l4);
  r3.subtract(v);
  affine S_inv (S);
  S_inv.invert_rigidbody();
  S_inv.mult(r3,r3_p);
  //S.print();S_inv.print();
  //cout<<"r4, r3, r3_p:"<<endl;r4.print();r3.print();r3_p.print();//exit(1);cout<<"s0="<<s0<<endl;
  // step 5
  limb_solver_yxx(r3_p.get_data(),q+1,ls+1,s0,bend);
  // step 6
  int ai[] = {2,1,0,0,0};
  chain_rotation(A,ai,q,4);
  //A.print();
  extvec ey (0,1,0), n2;
  A.mult(ey,n2);
  //n2.print();
  // step 7
  extvec ex (1,0,0), v1, v2;
  A.mult(ex,v1);
  n2.cross(n3,v2);
  //cout<<v1.dot(n3)<<endl;cout<<v1.dot(n2)<<endl;
  //cout<<asin(v1.dot(v2))<<endl;cout<<atan2(v1.dot(v2),n2.dot(n3))<<endl;
  q[4] = atan2(v1.dot(v2),n2.dot(n3))+(1-s0)*M_PI/2;
  //print_array(q,5);//exit(1);
  std::copy(q,q+5,jangles);
  //print_array(jangles,5);//exit(1);
  return true;
}


void redund_func0(int limbi, double* limb_rec);
void limbredund_func0(affine& A, const double* rec, double* constr);
void redund_gost_func0(int limbi, odepart* opart, double* limb_rec);

redundof::redundof(int index){
  dof = 3;
  redund_func = NULL;
  limbredund_func = NULL;
  redund_gost_func = NULL;
  switch (index) {
  case 4:
    dof = 4;//5;
    break;
  case 5:
    dof = 5;
    redund_func = &redund_func0;
    limbredund_func = &limbredund_func0;
    redund_gost_func = &redund_gost_func0;
    break;
  }
}

void redundof::set_rec_redundof(double* rec) const{
  if(redund_func == NULL){return;}
  double* p = rec+6;
  for(int i=0;i<n;i++){
    redund_func(i,p);
    p += dof;
  }
  //int rec_len = 6+dof*n;print_array(rec,rec_len);exit(1);
}

void redundof::constrlimb(affine& A, const double* rec, double* constr){
  if(limbredund_func == NULL){return;}
  limbredund_func(A,rec,constr);
}

void redundof::set_rec_redundof_gost(double* rec, vector<odepart*>& foot_oparts) const{
  if(redund_gost_func == NULL){return;}
  double* p = rec+6;
  for(int i=0;i<n;i++){
    redund_gost_func(i,foot_oparts[i],p);
    p += dof;
  }
}

// Extracts vector direction described by the spherical coordinate angles.
void vec_to_direction(const extvec& vec, double* direct){
  double x, y, z;
  vec.get_components(x,y,z);
  *direct = acos(z/vec.norm());
  *(direct+1) = atan2(y,x);
}

// Sets a unit vector along the direction direct
// described by two spherical coordinate angles.
void direction_to_vec(const double* direct, extvec& vec){
  double theta = *direct, phi = *(direct+1);
  double st = sin(theta);
  double x = st*cos(phi), y = st*sin(phi), z = cos(theta);
  vec.set(x,y,z);
  vec.set_v(3,0);
}

// Sets A to a rotation by angle angle around one of
// the coordinate axis axisi, where axisi = 0, 1 or 2
// for x, y or z axis respectively.
void set_xyz_rotation(affine& A, int axisi, double angle){
  double c = cos(angle), s = sin(angle);
  int i1 = (axisi+1)%3, i2 = (axisi+2)%3;
  A.set_unity();
  A.set_a(i1,i1,c);
  A.set_a(i2,i2,c);
  A.set_a(i1,i2,-s);
  A.set_a(i2,i1,s);
}

void chain_rotation(affine& A, int* axisi, double* q, int n){
  A.set_unity();
  affine B;
  for(int i=0;i<n;i++){
    set_xyz_rotation(B,axisi[i],q[i]);
    A.mult(B);
  }
}

void redund_func0(int limbi, double* limb_rec){
  double* p = limb_rec+3;
  for(int i=0;i<2;i++){*p++ = 0;}
  //*p++ = 0.2;*p = 1.57;
}

void limbredund_func0(affine& A, const double* rec, double* constr){
  extvec n, n1;
  direction_to_vec(rec,n);
  A.mult(n,n1);
  vec_to_direction(n1,constr);
  //n.print();n1.print();exit(1);
}

void redund_gost_func0(int limbi, odepart* opart, double* limb_rec){
  double* p = limb_rec+3;
  //for(int i=0;i<2;i++){*p++ = 0;}return;
  extvec foot_vec;
  opart->get_foot_vec_ground(foot_vec);
  foot_vec.times(-1);
  vec_to_direction(foot_vec,p);
  //print_array(p,2);foot_vec.print();
}

