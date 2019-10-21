#include "core.h"
#include "matrix.h"
#include "visualization.h"
#include "model.h"
#include "geom.h"

// computes rot that takes z to v; 
// rot is correct from the affine type perspective, 
// but needs to be transposed for the ode body rotation
void rot_ztov(dMatrix3& rot, extvec& v){
  extvec z (0,0,1), a;
  v.cross(z,a);
  double anorm = a.norm();
  if(anorm<1e-10){a.set(0,1,0);}
  double angle = asin(anorm/v.norm());
  if(v.dot(z)<0){angle = M_PI - angle;}
  double *p = a.get_data();
  dRFromAxisAndAngle(rot,*p,*(p+1),*(p+2),angle);

  // correctness check; leave for now, comment out later
  affine A; A.set_rotation(rot); extvec b,c; c.copy(v); c.normalize(); A.mult(z,b); c.subtract(b); if(c.norm()>1e-3){cout<<"A:"<<endl;A.print();cout<<"b:"<<endl;b.print();cout<<"v:"<<endl;v.print();cout<<"angle:"<<angle<<endl;cout<<"error: "<<b.norm()<<endl;exit(1);}
}

void transpose_odematrix(dMatrix3& m){
  dMatrix3 m1;
  dReal *p = m, *p1 = m1;
  for(int i=0;i<12;i++){*p1++ = *p++;}
  p = m;
  for(int i=0;i<3;i++){
    p1 = m1+i;
    for(int j=0;j<3;j++){
      *p++ = *p1;
      p1 += 4;
    }
    *p++ = 0;
  }
}

void posrot_from_affine(dVector3& pos, dMatrix3& rot, affine& A){
  double* a = A.get_data();
  dReal *p = pos, *r = rot;
  for(int i=0;i<12;i++){*r++ = *a++;}
  for(int i=0;i<3;i++){*p++ = *a++;}
}

void affine_from_posrot(affine& A, const dVector3& pos, const dMatrix3& rot){
  double* a = A.get_data();
  const dReal *p = pos, *r = rot;
  for(int i=0;i<12;i++){*a++ = *r++;}
  for(int i=0;i<3;i++){*a++ = *p++;}
  *a = 1;
}

void affine_from_orientation(affine& A, extvec* orientation){
  dVector3 pos;
  orientation[0].to_dvec(pos);
  dMatrix3 rot;
  double* p = orientation[1].get_data();
  dRFromEulerAngles(rot,*p,*(p+1),*(p+2));
  affine_from_posrot(A,pos,rot);
}

void mod_twopi(double& a){
  if (a < -M_PI) {
    while (a < -M_PI) {a += 2*M_PI;}
  } else if (a > M_PI) {
    while (a > M_PI) {a -= 2*M_PI;}
  }
}

// verify that euler angles are computed correctly
void euler_angles_from_affine(affine& A, double* angles){
  double r11, r21, r31, r32, r33;
  //double r12, r13;
  double *p = A.get_data();
  r11 = *p++; r21 = *p++; r31 = *p++; p++;
  p += 2; r32 = *p++; p++;
  p += 2; r33 = *p;

  double ps1, th1, ph1;
  th1 = -asin(r31);
  double ct1 = cos(th1);
  ps1 = atan2(r32/ct1,r33/ct1);
  ph1 = atan2(r21/ct1,r11/ct1);
  //cout <<"ps1="<<ps1<<" th1="<<th1<<" ph1="<<ph1<<endl;
  p = angles;
  *p++ = ps1;
  *p++ = th1;
  *p = ph1;

  if(fabs(ct1)<1e-5){cout<<"WARNING: small cos(theta)"<<endl;}
  return;

  /*double ps2, th2, ph2;
  th2 = M_PI - th1;
  double ct2 = cos(th2);
  ps2 = atan2(r32/ct2,r33/ct2);
  ph2 = atan2(r21/ct2,r11/ct2);
  //cout <<"ps2="<<ps2<<" th2="<<th2<<" ph2="<<ph2<<endl;*/
}

void print_dn(const dReal* a, int n){
  cout << "[";
  for(int i=0;i<4*n;i++){
    if (i%4) {cout << "\t";} else if (i) {cout << endl;}
    if (i) {cout << " ";}
    cout << a[i];
  }
  cout << "]" << endl;
}

void print_dmat(const dReal* mat){print_dn(mat,3);}

void print_dvec(const dReal* vec){print_dn(vec,1);}


visualizer::visualizer(){
  setup_odeworld();
  initialize_fn();
  player = NULL;
  view = new viewpoint;
  set_flag("manual_viewpoint",true);
  set_flag("texture",false);
  trimeshman = new trimeshmanager (odespace);
  //trimesh_test(); // temporary !!!!!!!!!!!1
}

visualizer::~visualizer(){
  unset_odeworld();
  delete view;
}

void visualizer::setup_odeworld(){
  dInitODE ();
  odeworld = dWorldCreate();
  odespace = dHashSpaceCreate(0);
  contact_group = dJointGroupCreate(0);
  
  double g = 1; //9.81;
  dWorldSetGravity (odeworld,0,0,-g);
  dWorldSetERP (odeworld,.8);
  //dWorldSetCFM (odeworld,1e-5);

  dCreatePlane (odespace,0,0,1,0);
  //dWorldSetQuickStepNumIterations(odeworld,100); // EXPERIMENTAL

}

void visualizer::unset_odeworld(){
  dJointGroupDestroy (contact_group);
  dSpaceDestroy (odespace);
  dWorldDestroy (odeworld);
  dCloseODE();
}

visualizer loop_vis;

void kinematicmodel::set_vis(){
  vis = &loop_vis;
  vis->model = this;
}

void vis_start(){
  //cout<<"starting visualizer"<<endl;
  loop_vis.set_viewpoint();
}

int speedup = 1;
void vis_loop(int){
  //cin.ignore();
  loop_vis.draw_inloop();
  for(int i=0;i<speedup;i++){loop_vis.step();}
  loop_vis.adjust_viewpoint();
}

/*
void vis_loop(int){
  //cin.ignore();
  loop_vis.draw_inloop();
  loop_vis.step();
  loop_vis.adjust_viewpoint();
}
*/
void visualizer::initialize_fn(){
  // setup pointers to drawstuff callback functions
  fn.version = DS_VERSION;
  fn.start = &vis_start;
  fn.step = &vis_loop;
  fn.stop = 0;
  fn.command = 0;
  fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;
}

void visualizer::start_loop(){
  char *argv1[3]; argv1[1] = new char[80];
  if(!texture_flag){strcpy(argv1[1],"-notex");}
  argv1[2] = new char[80];
  //strcpy(argv1[2],"-noshadow");
  //dsSimulationLoop (3,argv1,352,288,&fn);
  dsSimulationLoop (3,argv1,352*2,288*2,&fn);
}

void visualizer::draw(){
  //cout<<"starting loop"<<endl;
  start_loop();
}

void visualizer::draw_inloop(){
  int size = geoms.size();
  //cout<<size<<" geoms"<<endl;
  for(int i=0;i<size;i++){
    dGeomID geom = geoms[i];
    int geom_class = dGeomGetClass(geom);
    //cout<<"geomclass="<<geom_class<<endl;exit(1);
    dReal radius, length;

    switch (geom_class){
    case 0: {
      dReal rad = dGeomSphereGetRadius(geom);
      dsDrawSphere(dGeomGetPosition(geom),dGeomGetRotation(geom),rad);
    } break;
    case 1:
      dVector3 lengths;
      dGeomBoxGetLengths(geom,lengths);
      dsDrawBox(dGeomGetPosition(geom),dGeomGetRotation(geom),lengths);
      break;
    case 2:
      dGeomCapsuleGetParams(geom,&radius,&length);
      dsDrawCapsule(dGeomGetPosition(geom), dGeomGetRotation(geom), length, radius);
      break;
    case 3:
      dGeomCylinderGetParams(geom,&radius,&length);
      dsDrawCylinder(dGeomGetPosition(geom), dGeomGetRotation(geom), length, radius);
      break;
    case 8:
      trimeshman->draw(geom);
      break;
    default:
      cout << "unknown geom class " << geom_class << endl; exit(1);
      break;
    }
  }
  draw_forces();
}

void visualizer::set_viewpoint(){
  if(manual_viewpoint_flag){return;}
  float xyz[] = {0,0,0};
  float xyz_cam[] = {1,2,0.5};
  /*//float xyz_cam[] = {-1,-2,0.5};float hpr[] = {60,-20,0};
  float hpr[] = {245,-25,0};
  view->set(xyz,xyz_cam,hpr);*/
  view->set(xyz,xyz_cam);
}

void visualizer::adjust_viewpoint(){
  if(manual_viewpoint_flag){return;}
  const dReal* pos = dGeomGetPosition(geoms[0]);
  view->adjust(pos);
}

void visualizer::set_flag(string flag_name, bool value){
  if (flag_name == "manual_viewpoint") {
    manual_viewpoint_flag = value;
  } else if (flag_name == "texture") {
    texture_flag = value;
  } else if (flag_name == "smooth_view") {
    view->set_smooth(value);
  } else {
    cout << "ERROR: unknown flag " << flag_name << endl; exit(1);
  }
}

void nearCallback(void *data, dGeomID o1, dGeomID o2) {
  int cl[] = {dGeomGetClass(o1), dGeomGetClass(o2)};
  bool ignore_flag = true;
  for(int i=0;i<2;i++){
    if(cl[i]==dCapsuleClass || cl[i]==dSphereClass){
      if(cl[1-i]==dPlaneClass || cl[1-i]==dTriMeshClass){
	ignore_flag = false;
      }
    }
  }
  if(ignore_flag){return;}
  dBodyID b1 = dGeomGetBody(o1);
  dBodyID b2 = dGeomGetBody(o2);
  dContact contact;  
  contact.surface.mode = dContactBounce | dContactSoftCFM;
  // friction parameter
  contact.surface.mu = dInfinity; // experim
  //contact.surface.mu = 1;//10;//0;
  // bounce is the amount of "bouncyness".
  contact.surface.bounce = .5;//0.9;
  // bounce_vel is the minimum incoming velocity to cause a bounce
  contact.surface.bounce_vel = 0.1;
  // constraint force mixing parameter
  contact.surface.soft_cfm = 0.001;  
  //if (int numc = dCollide (o1,o2,1,&contact.geom,sizeof(dContact))) {
  if (dCollide (o1,o2,1,&contact.geom,sizeof(dContact))) {
    dJointID c = loop_vis.create_contact(&contact);
    dJointAttach (c,b1,b2);
  }
}

dJointID visualizer::create_contact(dContact* contact){
  return dJointCreateContact (odeworld, contact_group, contact);
}

void visualizer::simulate_odeworld(double dt_ode){
  dSpaceCollide(odespace,0,&nearCallback);
  dWorldQuickStep (odeworld, dt_ode);
  dJointGroupEmpty(contact_group);
}

void visualizer::set_speedup(int f){speedup = f;}

void visualizer::add_motor(dJointID hinge){
  motors.push_back(hinge);
}

void visualizer::set_ode_motor_torques(double* motor_torques){
  double *p = motor_torques;
  vector<dJointID>::iterator it = motors.begin();
  for(;it!=motors.end();it++){
    dJointAddHingeTorque(*it,*p++);
  }
}

void visualizer::get_ode_motor_angles(double* as){
  vector<dJointID>::iterator it = motors.begin();
  for(;it!=motors.end();it++){
    *as++ = dJointGetHingeAngle (*it);
  }
}

void visualizer::get_ode_motor_adas(double* as, double* das){
  vector<dJointID>::iterator it = motors.begin();
  for(;it!=motors.end();it++){
    dJointID hinge = (*it);
    *as++ = dJointGetHingeAngle (hinge);
    *das++ = dJointGetHingeAngleRate (hinge);
  }
}

void visualizer::get_ode_config(double* config){
  //odepart* torso_opart = (*model->get_odeparts())[0];
  odepart* torso_opart = get_torso_opart();
  affine A, B, C;
  torso_opart->get_frame_A_ground_from_body(A);
  modelnode* mnode = torso_opart->get_mnode();
  B.copy(*mnode->get_A_tobody());
  C.copy(*mnode->get_joint()->get_A_ground());
  B.invert_rigidbody();
  C.invert_rigidbody();
  C.mult(A);
  C.mult(B);
  //C.print();
  extvec pos;
  C.get_translation(pos);
  pos.get_components(config);
  euler_angles_from_affine(C,config+3);
  get_ode_motor_angles(config+6);
}

odepart* visualizer::get_torso_opart(){
  return (*model->get_odeparts())[0];
}

int draw_force_flag = 0;
void visualizer::add_force(dBodyID body, double* f){
  dBodyAddForce(body,f[0],f[1],f[2]);
  extvec force (f[0],f[1],f[2]);
  added_forces[body] = force;
  draw_force_flag = 5;
}

void visualizer::draw_forces(){
  switch(draw_force_flag){
  case 0: return; break;
  case 1: added_forces.clear(); break;
  default: draw_force_flag--; break;
  }
  map<dBodyID,extvec>::iterator it = added_forces.begin();
  for(;it!=added_forces.end();it++){
    dBodyID body = (*it).first;
    const dReal* pos0 = dBodyGetPosition(body);
    extvec pos (pos0[0],pos0[1],pos0[2]);
    extvec delpos = (*it).second;
    delpos.times(.01);
    pos.add(delpos);
    dVector3 pos1;
    pos.to_dvec(pos1);
    dsSetColor(0,1,0);
    dsDrawLine(pos0,pos1);
    dsDrawSphere(pos1,dBodyGetRotation(body),.1);
  }
}

void visualizer::trimesh_test(){
  for(int i=0;i<4;i++){
    dReal vert[3];
    for(int j=0;j<3;j++){vert[j] = rand()%10;}
    print_array<dReal>(vert,3);
    trimeshman->push_vertex(vert[0],vert[1],vert[2]);
  }
  trimeshman->push_triangle(0,1,2);
  trimeshman->push_triangle(0,1,3);
  dGeomID geom = trimeshman->new_trimesh();
  push_geom(geom);
}


void odepart::make(xml_node<>* xnode, modelnode* mnode_, visualizer* vis){
  mnode = mnode_;
  //part_name = xnode->first_attribute("name")->value();
  xmlnode_attr_to_val(xnode,"name",part_name);

  dWorldID world = *vis->get_odeworld();
  dSpaceID space = *vis->get_odespace();

  xml_node<>* geom_node = xnode->first_node("geom");
  string type = geom_node->first_attribute("type")->value();
  //cout<< type << endl;
  dBodyID body;
  if(type == "sphere"){
    //cout<<"making sphere"<<endl;
    double r, pos[3];
    xmlnode_attr_to_val(geom_node,"size",&r);
    xmlnode_attr_to_val(geom_node,"pos",pos);
    body = dBodyCreate(world);
    geom = dCreateSphere(space,r);
    dGeomSetBody(geom,body);
    A_geom.set_translation(pos);
    vis->push_geom(geom);
  } else if(type == "capsule"){
    make_ccylinder(vis,geom_node,true);
  } else if(type == "cylinder"){
    make_ccylinder(vis,geom_node,false);
  } else {
    cout << "WARNING: odepart of type " << type << " is undefined" << endl;
  }
}

void odepart::make_ccylinder(visualizer* vis, xml_node<>* geom_node, bool capped_flag){
  dWorldID world = *vis->get_odeworld();
  dSpaceID space = *vis->get_odespace();
  double r, fromto[6];
  xmlnode_attr_to_val(geom_node,"size",&r);
  xmlnode_attr_to_val(geom_node,"fromto",fromto);
  double len;
  dVector3 pos;
  dMatrix3 rot;
  capsule_lenposrot_from_fromto(len,pos,rot,fromto);
  dBodyID body = dBodyCreate(world);
  geom = (capped_flag)? dCreateCapsule(space,r,len) : dCreateCylinder(space,r,len);
  dGeomSetBody(geom,body);
  affine_from_posrot(A_geom,pos,rot);
  vis->push_geom(geom);
}

void odepart::capsule_lenposrot_from_fromto(double& len, dVector3& pos, dMatrix3& rot, double* fromto){
  for(int i=0;i<3;i++){pos[i]=(fromto[i]+fromto[i+3])/2.;}
  extvec r1, r2;
  r1.set(fromto);
  r2.set(fromto+3);
  r2.subtract(r1);
  len = r2.norm();
  rot_ztov(rot,r2);
  capsule_to_pos.set(fromto+3);
}

// computes ode body position and rotation from mnode body-frame
// note rot transposition to meet ODE convention
void odepart::get_body_posrot_from_frame(dVector3& pos, dMatrix3& rot){
  affine* A_frame = mnode->get_A_ground();
  affine A_body_ground;
  //cout<<"a_geom:"<<endl;A_geom.print_all();
  A_frame->mult(A_geom,A_body_ground);
  //cout<<"a_body_ground:"<<endl;A_body_ground.print_all();
  posrot_from_affine(pos,rot,A_body_ground);
  transpose_odematrix(rot);
}

void odepart::print(int detail_level){
  cout << "--- ode part ---" << endl;
  cout << "part name: " << part_name << endl;
  cout << "A_geom:" << endl;
  A_geom.print();

  if(detail_level > 1){print_ode();}
}

void odepart::print_ode(){
  dBodyID body = dGeomGetBody(geom);
  const dReal* pos = dBodyGetPosition(body);
  const dReal* vel = dBodyGetLinearVel(body);
  cout << "Position:" << endl;
  print_dvec(pos);
  cout << "LinearVel:" << endl;
  print_dvec(vel);
}

void odepart::get_com_pos(extvec& pos){
  extvec pos_body;
  A_geom.get_translation(pos_body);
  mnode->get_A_ground()->mult(pos_body,pos);
}

/*void odepart::get_foot_pos(extvec& pos){
  mnode->get_A_ground()->mult(capsule_to_pos,pos);
}
*/

void odepart::get_foot_pos(extvec& pos){
  get_foot_pos(pos, false);
}

void odepart::get_foot_pos(extvec& pos, bool from_body_flag){
  if (from_body_flag) {
    affine A;
    get_frame_A_ground_from_body(A);
    A.mult(capsule_to_pos, pos);
  } else {
    mnode->get_A_ground()->mult(capsule_to_pos,pos);
  } 
}

void odepart::make_fixed_joint(odepart* parent_part, visualizer* vis){
  dWorldID world = *vis->get_odeworld();
  dBodyID body = get_body();
  dBodyID parent_body = parent_part->get_body();
  dJointID joint = dJointCreateFixed(world,0);
  dJointAttach(joint,parent_body,body);
  dJointSetFixed(joint);
}

void odepart::make_hinge_joint(odepart* parent_part, visualizer* vis){
  dWorldID world = *vis->get_odeworld();
  dBodyID body = get_body();
  dBodyID parent_body = parent_part->get_body();

  modeljoint* joint = mnode->get_joint();
  affine* A = joint->get_A_ground();
  extvec pos, axis;
  A->get_translation(pos);
  axis.set(A->get_data() + 8);
  //cout<<"pos: ";pos.print();cout<<"axis: ";axis.print();cout<<endl;

  dJointID hinge = dJointCreateHinge(world,0);
  //dJointAttach(hinge,parent_body,body);
  dJointAttach(hinge,body,parent_body);
  double *u = pos.get_data();
  dJointSetHingeAnchor(hinge,u[0],u[1],u[2]);
  u = axis.get_data();
  dJointSetHingeAxis(hinge,u[0],u[1],u[2]);

  vis->add_motor(hinge);
}

// note rot transposition to meet ODE convention
void odepart::get_ode_body_A_ground(affine& A_ground){
  dVector3 pos;
  dMatrix3 rot;
  dBodyID body = get_body();
  const dReal* body_pos = dBodyGetPosition(body);
  const dReal* body_rot = dBodyGetRotation(body);
  dReal *p = pos, *p1 = rot;
  for(int i=0;i<4;i++){*p++ = *body_pos++;}
  for(int i=0;i<12;i++){*p1++ = *body_rot++;}
  transpose_odematrix(rot);
  affine_from_posrot(A_ground,pos,rot);
}

// computes A_ground of body-frame from ode body
void odepart::get_frame_A_ground_from_body(affine& A_ground){
  affine ode_A_ground;
  get_ode_body_A_ground(ode_A_ground);

  affine A_geom_inv;
  A_geom_inv.copy(A_geom);
  A_geom_inv.invert_rigidbody();
  ode_A_ground.mult(A_geom_inv,A_ground);

  // check: REMOVE LATER
  //affine A; A.copy(A_ground); A.subtract(*mnode->get_A_ground()); float norm = A.norm(); if(norm>1e-5){cout << "ERROR: norm = " << norm << endl; exit(1);}
}


viewpoint::viewpoint(){
  k0 = .0002;
  for(int i=0;i<3;i++){
    xyz0[i] = 0;
    xyz_rate[i] = 0;
    xyz_ref_rate[i] = 0;    
  }
  smooth_flag = false;
}

void viewpoint::set(float* xyz, float* xyz_cam, float* hpr_){
  set_xyz(xyz, xyz_cam);
  for(int i=0;i<3;i++){hpr[i] = hpr_[i];}
}

void viewpoint::set(float* xyz, float* xyz_cam){
  set_xyz(xyz, xyz_cam);
  hpr_from_cam_rel();
}

void viewpoint::adjust(const dReal* xyz){
  if (smooth_flag) {smooth_xyzref_update(xyz);}
  else {hard_xyzref_update(xyz);}
  float xyz_cam[3];
  for(int i=0;i<3;i++){
    xyz_cam[i] = xyz_ref[i] + xyz_cam_rel[i];
  }
  dsSetViewpoint(xyz_cam,hpr);
}

void viewpoint::hard_xyzref_update(const dReal* xyz){
  for(int i=0;i<3;i++){xyz_ref[i] = xyz[i];}
}

void viewpoint::smooth_xyzref_update(const dReal* xyz){
  float k = k0*speedup;
  for(int i=0;i<3;i++){
    xyz_ref[i] += xyz_ref_rate[i];
    xyz_ref_rate[i] -= k*(xyz_ref[i]-xyz[i])+2*sqrt(k)*(xyz_ref_rate[i]-xyz_rate[i]);
    xyz_rate[i] = xyz[i] - xyz0[i];
    xyz0[i] = xyz[i];
  }
}

void viewpoint::print(){
  cout << "----- veiwpoint -----" << endl;
  print_array(xyz_ref,3,"xyz_ref: ");
  print_array(xyz_cam_rel,3,"xyz_cam_rel: ");
  print_array(hpr,3,"hpr: ");
}

void viewpoint::hpr_from_cam_rel(){
  float *p = xyz_cam_rel, *p1 = hpr;
  float x, y, z;
  x = *p++; y = *p++; z = *p;
  float l = sqrt(x*x+y*y);
  float cf = 180/M_PI; // cf = conversion factor
  *p1++ = atan2(-y,-x)*cf;
  *p1++ = asin((-z-.45)/l)*cf;
  *p1 = 0;
}

void viewpoint::set_xyz(float* xyz, float* xyz_cam){
  for(int i=0;i<3;i++){
    xyz_ref[i] = xyz[i];
    xyz_cam_rel[i] += xyz_cam[i];
  }
}

void viewpoint::shift_cam(float x, float y, float z){
  float del_xyz_cam[] = {x, y, z};
  for(int i=0;i<3;i++){
    xyz_cam_rel[i] += del_xyz_cam[i];
  }
}
