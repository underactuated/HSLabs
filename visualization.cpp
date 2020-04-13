#include "visualization.h"
#include "model.h"
#include "geom.h"

//TODO: support for more geom classes

// Computes rot that takes z to v; 
// rot is correct from the affine type perspective, 
// but needs to be transposed for the ode-body rotation.
// ode-body = ODE body
void rot_ztov(dMatrix3& rot, const extvec& v){
  //extvec z (0,0,1), a;
  const extvec z (0,0,1);
  extvec a;
  v.cross(z,a);
  double anorm = a.norm();
  if(anorm<1e-10){a.set(0,1,0);}
  double angle = asin(anorm/v.norm());
  if(v.dot(z)<0){angle = M_PI - angle;}
  double *p = a.get_data();
  dRFromAxisAndAngle(rot,*p,*(p+1),*(p+2),angle);

  // correctness check; TODO: leave for now, comment out later
  affine A; A.set_rotation(rot); extvec b,c; c=v; c.normalize(); A.mult(z,b); c.subtract(b); if(c.norm()>1e-3){cout<<"A:"<<endl;A.print();cout<<"b:"<<endl;b.print();cout<<"v:"<<endl;v.print();cout<<"angle:"<<angle<<endl;cout<<"error: "<<b.norm()<<endl;exit(1);}
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

// Initializes translation pos and rotation rot 
// from an affine transformation A. Note: what we 
// call an affine transformation in this project,
// is really a rigid body transformation. 
void posrot_from_affine(dVector3& pos, dMatrix3& rot, const affine& A){
  const double* a = A.get_data();
  dReal *p = pos, *r = rot;
  for(int i=0;i<12;i++){*r++ = *a++;}
  for(int i=0;i<3;i++){*p++ = *a++;}
}

// reverse of posrot_from_affine
void affine_from_posrot(affine& A, const dVector3& pos, const dMatrix3& rot){
  double* a = A.get_data();
  const dReal *p = pos, *r = rot;
  for(int i=0;i<12;i++){*a++ = *r++;}
  for(int i=0;i<3;i++){*a++ = *p++;}
  *a = 1;
}

void affine_from_orientation(affine& A, const extvec* orientation){
  dVector3 pos;
  orientation[0].to_dvec(pos);
  dMatrix3 rot;
  const double* p = orientation[1].get_data();
  dRFromEulerAngles(rot,*p,*(p+1),*(p+2));
  affine_from_posrot(A,pos,rot);
}

// analog of arrayops::modulus(a,2*M_PI)
// that works for any a
void mod_twopi(double& a){
  if (a < -M_PI) {
    while (a < -M_PI) {a += 2*M_PI;}
  } else if (a > M_PI) {
    while (a > M_PI) {a -= 2*M_PI;}
  }
}

void euler_angles_from_affine(const affine& A, double* angles){
  double r11, r21, r31, r32, r33;
  const double *p = A.get_data();
  r11 = *p++; r21 = *p++; r31 = *p++; p++;
  p += 2; r32 = *p++; p++;
  p += 2; r33 = *p;

  double ps1, th1, ph1;
  th1 = -asin(r31);
  double ct1 = cos(th1);
  ps1 = atan2(r32/ct1,r33/ct1);
  ph1 = atan2(r21/ct1,r11/ct1);
  //cout <<"ps1="<<ps1<<" th1="<<th1<<" ph1="<<ph1<<endl;
  double *p1 = angles;
  *p1++ = ps1;
  *p1++ = th1;
  *p1 = ph1;

  if(fabs(ct1)<1e-5){cout<<"WARNING: small cos(theta)"<<endl;}
  return;
}

// for print_dmat and print_dvec
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
  //trimesh_test(); // temporary !!!
  speedup = 1;
}

visualizer::~visualizer(){
  unset_odeworld();
  delete view;
  delete trimeshman;
}

// note the gravity constant g = 1
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

// loop_vis is used by various free functions defined in this file
// that are used by ODE and drawstuff callback functions, e.g.
// vis_start(), vis_loop(), nearCallback().
// See also visualizer::initialize_fn
static visualizer loop_vis; 

// Sets loop_vis as model's visualizer.
void kinematicmodel::set_vis(){
  vis = &loop_vis;
  vis->set_model(this);
}

// used by fn
void vis_start(){
  loop_vis.set_viewpoint();
}

// used by fn
void vis_loop(int){
  //cin.ignore();
  int speedup = loop_vis.get_speedup();
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

// Initializes drawstuff (ds) callback function fn.
void visualizer::initialize_fn(){
  // setup pointers to drawstuff callback functions
  fn.version = DS_VERSION;
  fn.start = &vis_start;
  fn.step = &vis_loop;
  fn.stop = 0;
  fn.command = 0;
  fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;
}

// Starts ds simulation loop, with possible no-texture 
// and/or no-shadow flags; specifies ds window size.
void visualizer::start_loop(){
  char *argv1[3]; argv1[1] = new char[80];
  if(!texture_flag){strcpy(argv1[1],"-notex");}
  argv1[2] = new char[80];
  //strcpy(argv1[2],"-noshadow");
  //dsSimulationLoop (3,argv1,352,288,&fn);
  dsSimulationLoop (3,argv1,352*2,288*2,&fn);
}

// Launches ds simulation/drawing loop.
void visualizer::draw(){
  //cout<<"starting loop"<<endl;
  start_loop();
}

// Draws all geoms in every iteration of ds loop.
// Current implementation includes:
// sphere, box, capsule, cylinder, trimesh
void visualizer::draw_inloop(){
  int size = geoms.size();
  for(int i=0;i<size;i++){
    dGeomID geom = geoms[i];
    int geom_class = dGeomGetClass(geom);
    //cout<<"geomclass="<<geom_class<<endl;exit(1);
    dReal radius, length;

    switch (geom_class){
    case dSphereClass: {
      dReal rad = dGeomSphereGetRadius(geom);
      dsDrawSphere(dGeomGetPosition(geom),dGeomGetRotation(geom),rad);
    } break;
    case dBoxClass:
      dVector3 lengths;
      dGeomBoxGetLengths(geom,lengths);
      dsDrawBox(dGeomGetPosition(geom),dGeomGetRotation(geom),lengths);
      break;
    case dCapsuleClass:
      dGeomCapsuleGetParams(geom,&radius,&length);
      dsDrawCapsule(dGeomGetPosition(geom), dGeomGetRotation(geom), length, radius);
      break;
    case dCylinderClass:
      dGeomCylinderGetParams(geom,&radius,&length);
      dsDrawCylinder(dGeomGetPosition(geom), dGeomGetRotation(geom), length, radius);
      break;
    case dTriMeshClass:
      trimeshman->draw(geom);
      break;
    default:
      cout << "unknown geom class " << geom_class << endl; exit(1);
      break;
    }
  }
  draw_forces();
}

// Sets initial camera's viewpoint.
void visualizer::set_viewpoint(){
  if(manual_viewpoint_flag){return;}
  float xyz[] = {0,0,0};
  float xyz_cam[] = {1,2,0.5};
  /*//float xyz_cam[] = {-1,-2,0.5};float hpr[] = {60,-20,0};
  float hpr[] = {245,-25,0};
  view->set(xyz,xyz_cam,hpr);*/
  view->set(xyz,xyz_cam);
}

// Adjusts camera's viewpoint on every ds loop iteration.
void visualizer::adjust_viewpoint(){
  if(manual_viewpoint_flag){return;}
  const dReal* pos = dGeomGetPosition(geoms[0]);
  view->adjust(pos);
}

void visualizer::set_flag(string flag_name, bool value){
  if (flag_name == "manual_viewpoint") {  // for manual adjustment of camera
    manual_viewpoint_flag = value;
  } else if (flag_name == "texture") {  // for showing sky and ground texture
    texture_flag = value;
  } else if (flag_name == "smooth_view") {
    // if smooth_view is not set, camera is locked on torso com
    view->set_smooth(value);
  } else {
    cout << "ERROR: unknown flag " << flag_name << endl; exit(1);
  }
}

// a callback function for computing geom collisions
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
  contact.surface.bounce = 0.5;//0.9;
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

// Creates ODE contact on collision.
dJointID visualizer::create_contact(dContact* contact){
  return dJointCreateContact (odeworld, contact_group, contact);
}

// simulation step of length dt_ode
void visualizer::simulate_odeworld(double dt_ode){
  dSpaceCollide(odespace,0,&nearCallback);
  dWorldQuickStep (odeworld, dt_ode);
  dJointGroupEmpty(contact_group);
}

// Sets a speedup factor relative to nominal ds visualization rate.
void visualizer::set_speedup(int f){
  speedup = f;
  view->set_speedup(speedup);
}

// Stores a motorized joint.
void visualizer::add_motor(dJointID hinge){
  motors.push_back(hinge);
}

void visualizer::set_ode_motor_torques(const double* motor_torques){
  const double *p = motor_torques;
  vector<dJointID>::iterator it = motors.begin();
  for(;it!=motors.end();it++){
    dJointAddHingeTorque(*it,*p++);
  }
}

// Gets motorized joint coordinates (angles).
void visualizer::get_ode_motor_angles(double* as) const {
  vector<dJointID>::const_iterator it = motors.begin();
  for(;it!=motors.end();it++){
    *as++ = dJointGetHingeAngle (*it);
  }
}

// Gets motor's angles and angle rates (as and das respectively).
void visualizer::get_ode_motor_adas(double* as, double* das) const {
  vector<dJointID>::const_iterator it = motors.begin();
  for(;it!=motors.end();it++){
    dJointID hinge = (*it);
    *as++ = dJointGetHingeAngle (hinge);
    *das++ = dJointGetHingeAngleRate (hinge);
  }
}

// Computes model's configuration config from 
// torso ode part's location and motor angles.
void visualizer::get_ode_config(double* config) const {
  const odepart* torso_opart = get_torso_opart();
  affine A;
  //torso_opart->get_frame_A_ground_from_body(A);
  torso_opart->get_A_ground_body_from_odebody(A);
  modelnode* mnode = torso_opart->get_mnode();
  affine B (*mnode->get_A_pj_body());
  affine C (*mnode->get_joint()->get_A_ground());
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

// opart = ode part
const odepart* visualizer::get_torso_opart() const {
  return model->get_odepart(0);
}

int draw_force_flag = 0;
// Adds a perturbing force to ode-body 
// and stores it for ds visualization.
void visualizer::add_force(dBodyID odebody, const double* f){
  dBodyAddForce(odebody,f[0],f[1],f[2]);
  extvec force (f[0],f[1],f[2]);
  added_forces[odebody] = force;
  draw_force_flag = 5; // number of ds frames that visualization persists
}

// Draws perturbing forces.
// TODO: replace sphere with cone for arrow tip
void visualizer::draw_forces(){
  switch(draw_force_flag){
  case 0: return; break;
  case 1: added_forces.clear(); break;
  default: draw_force_flag--; break;
  }
  map<dBodyID,extvec>::iterator it = added_forces.begin();
  for(;it!=added_forces.end();it++){
    dBodyID odebody = (*it).first;
    const dReal* pos0 = dBodyGetPosition(odebody);
    extvec pos (pos0[0],pos0[1],pos0[2]);
    extvec delpos = (*it).second;
    delpos.times(.01);
    pos.add(delpos);
    dVector3 pos1;
    pos.to_dvec(pos1);
    dsSetColor(0,1,0);
    dsDrawLine(pos0,pos1);
    dsDrawSphere(pos1,dBodyGetRotation(odebody),.1);
  }
}


// Creates an ODE part (body and geom) corresponding to 
// xml_node/modelnode (xnode and mnode, respectively).
// Geom is stored in visualizer vis.
// Currently supported geoms: sphere, capsule, cylinder
void odepart::make(const xml_node<>* xnode, modelnode* mnode_, visualizer* vis){
  mnode = mnode_;
  //part_name = xnode->first_attribute("name")->value();
  xmlnode_attr_to_val(xnode,"name",part_name);

  dWorldID world = *vis->get_odeworld();
  dSpaceID space = *vis->get_odespace();

  xml_node<>* geom_node = xnode->first_node("geom");
  string type = geom_node->first_attribute("type")->value();
  //cout<< type << endl;
  dBodyID odebody;
  if(type == "sphere"){
    double r, pos[3];
    xmlnode_attr_to_val(geom_node,"size",&r);
    xmlnode_attr_to_val(geom_node,"pos",pos);
    odebody = dBodyCreate(world);
    geom = dCreateSphere(space,r);
    dGeomSetBody(geom,odebody);
    A_body_geom.set_translation(pos);
    vis->push_geom(geom);
  } else if(type == "capsule"){
    make_ccylinder(vis,geom_node,true);
  } else if(type == "cylinder"){
    make_ccylinder(vis,geom_node,false);
  } else {
    cout << "WARNING: odepart of type " << type << " is undefined" << endl;
  }
}

// Makes capsule or cylinder.
// If capped_flag, makes capsule.
// Used by make().
void odepart::make_ccylinder(visualizer* vis, const xml_node<>* geom_node, bool capped_flag){
  dWorldID world = *vis->get_odeworld();
  dSpaceID space = *vis->get_odespace();
  double r, fromto[6];
  xmlnode_attr_to_val(geom_node,"size",&r);
  xmlnode_attr_to_val(geom_node,"fromto",fromto);
  double len;
  dVector3 pos;
  dMatrix3 rot;
  capsule_lenposrot_from_fromto(len,pos,rot,fromto);
  dBodyID odebody = dBodyCreate(world);
  geom = (capped_flag)? dCreateCapsule(space,r,len) : dCreateCylinder(space,r,len);
  dGeomSetBody(geom,odebody);
  affine_from_posrot(A_body_geom,pos,rot);
  vis->push_geom(geom);
  if(capped_flag){rcap = r;}
}

// Computes capsule's length, position and rotation (in body frame)
// from fromto parameters given in xml file.
void odepart::capsule_lenposrot_from_fromto(double& len, dVector3& pos, dMatrix3& rot, const double* fromto){
  for(int i=0;i<3;i++){pos[i]=(fromto[i]+fromto[i+3])/2.;}
  extvec r1, r2;
  r1.set(fromto);
  r2.set(fromto+3);
  r2.subtract(r1);
  len = r2.norm();
  rot_ztov(rot,r2);
  capsule_to_pos.set(fromto+3);
}

// Computes ode-body position and rotation from mnode body-frame.
// Note rot transposition to meet ODE convention.
//void odepart::get_body_posrot_from_frame(dVector3& pos, dMatrix3& rot){
void odepart::get_odebody_posrot_from_body(dVector3& pos, dMatrix3& rot){
  const affine* A_ground_body = mnode->get_A_ground();
  affine A_ground_odebody;
  //cout<<"a_geom:"<<endl;A_body_geom.print_all();
  A_ground_body->mult(A_body_geom,A_ground_odebody);
  //cout<<"A_ground_odebody:"<<endl;A_ground_odebody.print_all();
  posrot_from_affine(pos,rot,A_ground_odebody);
  transpose_odematrix(rot);
}

void odepart::print(int detail_level){
  cout << "--- ode part ---" << endl;
  cout << "part name: " << part_name << endl;
  cout << "A_body_geom:" << endl;
  A_body_geom.print();

  if(detail_level > 1){print_ode();}
  if(detail_level > 2){mnode->print(detail_level-3);}
}

// Prints ode-body position and velocity.
void odepart::print_ode(){
  dBodyID odebody = dGeomGetBody(geom);
  const dReal* pos = dBodyGetPosition(odebody);
  const dReal* vel = dBodyGetLinearVel(odebody);
  cout << "Position:" << endl;
  print_dvec(pos);
  cout << "LinearVel:" << endl;
  print_dvec(vel);
}

// Computes com from model node.
void odepart::get_com_pos(extvec& pos) const {
  extvec pos_body;
  A_body_geom.get_translation(pos_body);
  mnode->get_A_ground()->mult(pos_body,pos);
}

/*void odepart::get_foot_pos(extvec& pos){
  mnode->get_A_ground()->mult(capsule_to_pos,pos);
}
*/

// Computes foot position from model node.
void odepart::get_foot_pos(extvec& pos) const {
  get_foot_pos(pos, false);
}

// Computes foot position from ode-body if from_body_flag,
// otherwise from model node.
void odepart::get_foot_pos(extvec& pos, bool from_body_flag) const {
  if (from_body_flag) {
    affine A;
    //get_frame_A_ground_from_body(A);
    get_A_ground_body_from_odebody(A);
    A.mult(capsule_to_pos, pos);
  } else {
    mnode->get_A_ground()->mult(capsule_to_pos,pos);
  } 
}

// Makes a fixed ODE joint (when no joint is specified 
// in xml file for a given body).
void odepart::make_fixed_joint(odepart* parent_part, visualizer* vis){
  dWorldID world = *vis->get_odeworld();
  dBodyID odebody = get_odebody();
  dBodyID parent_odebody = parent_part->get_odebody();
  dJointID joint = dJointCreateFixed(world,0);
  dJointAttach(joint,parent_odebody,odebody);
  dJointSetFixed(joint);
}

// Makes a hinge ODE joint, stores it visualizer.
// For now any hinge is assumed motorized.
void odepart::make_hinge_joint(odepart* parent_part, visualizer* vis){
  dWorldID world = *vis->get_odeworld();
  dBodyID odebody = get_odebody();
  dBodyID parent_odebody = parent_part->get_odebody();

  modeljoint* joint = mnode->get_joint();
  affine* A = joint->get_A_ground();
  extvec pos, axis;
  A->get_translation(pos);
  axis.set(A->get_data() + 8);
  //cout<<"pos: ";pos.print();cout<<"axis: ";axis.print();cout<<endl;

  dJointID hinge = dJointCreateHinge(world,0);
  dJointAttach(hinge,odebody,parent_odebody);
  double *u = pos.get_data();
  dJointSetHingeAnchor(hinge,u[0],u[1],u[2]);
  u = axis.get_data();
  dJointSetHingeAxis(hinge,u[0],u[1],u[2]);

  vis->add_motor(hinge);
}

// Computes ode-body/geom transformation A_ground relative to ground frame.
// Note rot transposition to meet ODE convention.
//void odepart::get_ode_body_A_ground(affine& A_ground) const {
void odepart::get_A_ground_odebody(affine& A_ground) const {
  dVector3 pos;
  dMatrix3 rot;
  dBodyID odebody = get_odebody();
  const dReal* odebody_pos = dBodyGetPosition(odebody);
  const dReal* odebody_rot = dBodyGetRotation(odebody);
  dReal *p = pos, *p1 = rot;
  for(int i=0;i<4;i++){*p++ = *odebody_pos++;}
  for(int i=0;i<12;i++){*p1++ = *odebody_rot++;}
  transpose_odematrix(rot);
  affine_from_posrot(A_ground,pos,rot);
}

// Computes A_ground of body-frame from ode-body.
// TODO: distinction between body and ode-body (hence between body-frame, and ode-body-frame) needs to be clarified.
//void odepart::get_frame_A_ground_from_body(affine& A_ground) const {
void odepart::get_A_ground_body_from_odebody(affine& A_ground) const {
  affine ode_A_ground;
  get_A_ground_odebody(ode_A_ground);

  affine A_body_geom_inv (A_body_geom);
  A_body_geom_inv.invert_rigidbody();
  ode_A_ground.mult(A_body_geom_inv,A_ground);

  // check: REMOVE LATER
  //affine A; A=A_ground; A.subtract(*mnode->get_A_ground()); float norm = A.norm(); if(norm>1e-5){cout << "ERROR: norm = " << norm << endl; exit(1);}
}

// Computes foot vector (capsule_to_pos with ext component = 0)
// in the ground frame, using foot ODE body.
void odepart::get_foot_vec_ground(extvec& vec) const {
  affine A;
  get_A_ground_body_from_odebody(A);
  extvec vec0 (capsule_to_pos);
  vec0.set_v(3,0);
  A.mult(vec0,vec);
}


// Camera position is either adjusted manually if manual_viewpoint_flag
// (then viewpoint plays no role), or it is adjusted automatically
// on every ds loop iteration. If smooth_flag the automatic adjustment
// is smooth, otherwise it is hard. In smooth adjustment the camera
// is positioned relative to reference point, that smoothly tracks
// torso com (via a DP control). In hard adjustment reference point
// coinsides with torso com.
// xyz and xyz_rate are torso position and velocity.
// xyz_ref and xyz_ref_rate are reference point and velocity, relative
// to which camera's position is defined.
// k0 is the gain of the tracking (critical) PD controller.
viewpoint::viewpoint(){
  k0 = .0002;
  for(int i=0;i<3;i++){
    xyz0[i] = 0;
    xyz_rate[i] = 0;
    xyz_ref_rate[i] = 0;    
  }
  smooth_flag = false;
  speedup = 1;
  scale = 1;
}

// Sets initial reference point xyz and camera position relative
// to reference point xyz_cam, and camera orientation
// hpr = heading, pitch, roll. 
void viewpoint::set(const float* xyz, const float* xyz_cam, const float* hpr_){
  set_xyz(xyz, xyz_cam);
  for(int i=0;i<3;i++){hpr[i] = hpr_[i];}
}

// Sets initial reference point xyz and camera position relative
// to reference point xyz_cam, camera orientation is computed
// so it faces reference point.
void viewpoint::set(const float* xyz, const float* xyz_cam){
  set_xyz(xyz, xyz_cam);
  hpr_from_cam_rel();
}

// Adjusts camera position, taking com position as argument.
void viewpoint::adjust(const dReal* xyz){
  if (smooth_flag) {smooth_xyzref_update(xyz);}
  else {hard_xyzref_update(xyz);}
  float xyz_cam[3];
  for(int i=0;i<3;i++){
    //xyz_cam[i] = xyz_ref[i] + xyz_cam_rel[i];
    xyz_cam[i] = xyz_ref[i] + xyz_cam_rel[i]*scale;
  }
  dsSetViewpoint(xyz_cam,hpr);
}

// Updates xyz_ref.
void viewpoint::hard_xyzref_update(const dReal* xyz){
  for(int i=0;i<3;i++){xyz_ref[i] = xyz[i];}
}

void viewpoint::smooth_xyzref_update(const dReal* xyz){
  double k = k0*speedup;
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

// Computes camera orientation from its relative position.
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

// Sets initial reference point and relative camera position.
void viewpoint::set_xyz(const float* xyz, const float* xyz_cam){
  for(int i=0;i<3;i++){
    xyz_ref[i] = xyz[i];
    xyz_cam_rel[i] += xyz_cam[i];
  }
}

// Shifts relative camera position.
void viewpoint::shift_cam(float x, float y, float z){
  float del_xyz_cam[] = {x, y, z};
  for(int i=0;i<3;i++){
    xyz_cam_rel[i] += del_xyz_cam[i];
  }
}
