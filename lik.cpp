#include "lik.h"


// Initializes LIK solver for a given model.
// Solver functions must be implemented in 
// this file for every model.
liksolver::liksolver(const kinematicmodel* model){
  string xmlfname = model->get_xmlfname();
  if(xmlfname == "myant.xml"){index = 0;}
  else if(xmlfname == "hexapod.xml"){index = 1;}
  else if(xmlfname == "spider.xml"){index = 2;}
  else {
    index = -1;
    cout << "WARNING: limb inverse kinematics (LIK) solver is not set for " << xmlfname << endl;
    return;
  }
  solver_func = NULL;
  rcap = 0;
  set_limbs(model);
}

liksolver::~liksolver(){
  while(limbs.size()){
    delete limbs.back();
    limbs.pop_back();
  }
  if(solver_func){delete [] solver_func;}
}

// limb_solver# solves for joint_angles, given limbi, pos_limb and (knee) bend
// returns true if solution exists.
bool limb_solver0(int limbi, extvec& pos_limb, extvec& joint_angles, bool bend);
bool limb_solver1(int limbi, extvec& pos_limb, extvec& joint_angles, bool bend);
bool limb_solver2(int limbi, extvec& pos_limb, extvec& joint_angles, bool bend);

// bend_solver# solves for and returns bend, given limbi and angles
// optionally solving for (foot) pos, if pos_flag set true.
bool bend_solver0(int limbi, extvec& pos, extvec& angles, bool pos_flag);
bool bend_solver1(int limbi, extvec& pos, extvec& angles, bool pos_flag);
bool bend_solver2(int limbi, extvec& pos, extvec& angles, bool pos_flag);

// Sets up liklimb objects (limbs) containing a lik solver for every limb.
// Array inds lists model node indices of top-link limb nodes.
void liksolver::set_limbs(const kinematicmodel* model){
  //model->print();
  vector<int> child_inds;
  solver_func = new SolverFuncType [2];
  switch (index) {
  case 0: {
    int inds[] = {2,6,10,14};
    child_inds.insert(child_inds.end(),inds,inds+4);
    solver_func[0] = &limb_solver0;
    solver_func[1] = &bend_solver0;
  } break;
  case 1: {
    int inds[] = {2,5,9,12,16,19};
    child_inds.insert(child_inds.end(),inds,inds+6);
    solver_func[0] = &limb_solver1;
    solver_func[1] = &bend_solver1;
  } break;
  case 2: {
    int inds[] = {1,4,7,10,13,16};
    child_inds.insert(child_inds.end(),inds,inds+6);
    solver_func[0] = &limb_solver2;
    solver_func[1] = &bend_solver2;
  } break;
  default:
    cout << "WARNING: case " << index << " undefined in set_limbs" << endl;
    break;
  }
  int imax = child_inds.size();
  for(int i=0;i<imax;i++){
    const modelnode* child = model->get_mnode(child_inds[i]);
    limbs.push_back(new liklimb (i,child));
    limbs.back()->set_solver_func(solver_func);
  }
  set_rcap(model,child_inds);
}

// Arranges a limb (indexed by limbi) in the model (by setting joint values)
// so that its foot is placed at (x,y,z) (in absolute coordinates). 
void liksolver::place_limb(int limbi, double x, double y, double z) const {
  extvec pos_ground (x,y,z);
  limbs[limbi]->place_limb(pos_ground);
}

// rec contains torso orientation (first 6 components)
// and ground positions for all feet (the rest 3*imax components).
void liksolver::place_limbs(const double* rec) const {
  if(index == -1){cout << "WARNING: LIK undefined" << endl; return;}
  const double* p = rec;
  int imax = limbs.size();
  for(int i=0;i<imax;i++){
    extvec pos_ground;
    pos_ground.set(p);
    limbs[i]->place_limb(pos_ground);
    p += 3;
  }
}

// Gets limb position (body position of top-link model node),
// that usually coinsides with the hip-joint position.
// Therefore, we also referr to it as hip-position.
void liksolver::get_limb_pos0(int limbi, extvec& pos) const {
  limbs[limbi]->get_pos0(pos);
}

void liksolver::print_limb_pos0s() const {
  cout << "--- limb hip positions ---" << endl;
  for(int i=0;i<(int)limbs.size();i++){
    extvec pos0;
    get_limb_pos0(i,pos0);
    cout << i << ": ";
    pos0.print();
  }
}

// Tests LIK solver functions. Currently only tests
// consistency of lik_solver and bend_solver funcsions,
// and only for yxx join axes configurations.
// TODO: implement testing for other limb configurations,
// verification of lik_solver correctness using model.
void liksolver::solver_test(int n) const {
  cout << "LIK solver testing ..." << endl;
  vector<liklimb*>::const_iterator it = limbs.begin();
  for(;it!=limbs.end();it++){(*it)->solver_test_yxx(n);}
  cout << "... success!" << endl;
}

// Extracts rcap from foot odeparts, and verifies that
// all feet have capsules of the same size.
void liksolver::set_rcap(const kinematicmodel* model, const vector<int>& limb_inds){
  if(!model->if_vis()){return;}
  vector<int>::const_iterator it = limb_inds.begin();
  for(;it!=limb_inds.end();it++){
    double rcap1 = model->get_odepart((*it)+2)->get_rcap();
    if((rcap > 0) && (rcap1 != rcap)){cout<<"ERROR: currently, rcaps must be same for all feet, rcap = "<<rcap<<", rcap1 = "<<rcap1<<endl;exit(1);}
    rcap = rcap1;
  }
}


// Solves 3-link limb with hinge axis y-x-x, (quadruped/hexapod configuration).
// Foot position pos_limb is in the hip joint's frame.
// To change knee bend (inward/outward) flip sign of s1.
bool limb_solver_yxx(int limbi, const extvec& pos_limb, extvec& joint_angles, const double* ls, int ysign, bool bend){
  double l0=ls[0], l1=ls[1], l2=ls[2];

  int s0 = ysign;
  int s1 = 2*int(bend)-1;

  extvec pos0 (0,0,s0*l0);
  extvec pos1 (pos_limb);
  pos1.subtract(pos0);
  double l = pos1.norm();
  if(l1+l2-l<0){cout<<"LIK ERROR: limb position is unreachable"<<endl;return false;}

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

  //joint_angles.print();
  return true;
}

// Solves 3-link limb with hinge axis z-x-x, (spider configuration).
// Foot position pos_limb is in the hip joint's frame.
// To change knee bend (inward/outward) flip sign of s1.
bool limb_solver_zxx(int limbi, const extvec& pos_limb, extvec& joint_angles, const double* ls, int ysign, bool bend){
  double l0=ls[0], l1=ls[1], l2=ls[2];

  int s0 = ysign;
  int s1 = 2*int(bend)-1;

  extvec pos0 (0,0,l0);
  extvec pos1 (pos_limb);
  pos1.add(pos0);
  double l = pos1.norm();
  if(l1+l2-l<0){cout<<"LIK ERROR: limb position is unreachable"<<endl;return false;}

  double x, y, z;
  pos_limb.get_components(x,y,z);
  double c = (z+l0)/l;
  double phi = atan2(x,y);
  //double theta = acos(c) + (1-s0)*M_PI/2;
  double theta = acos(c) - s0*M_PI/2;

  mod_twopi(phi);
  mod_twopi(theta);

  double ll = l*l;
  double del = l2*l2-l1*l1;
  double beta = s1*acos((ll-del)/(2*l1*l));
  double gamma = s1*acos((ll+del)/(2*l2*l));
  joint_angles.set(-phi,-theta+beta,-(beta+gamma));

  //joint_angles.print();
  return true;
}

// limb link lengths, (hex/quad and spider configurations).
double limb_ls[] = {.05, .4, .4};
double limb_ls1[] = {.1, .4, .4};

// quadruped lik solver
bool limb_solver0(int limbi, extvec& pos_limb, extvec& joint_angles, bool bend){
  int ysign = (limbi<2)? 1 : -1;
  return limb_solver_yxx(limbi,pos_limb,joint_angles,limb_ls,ysign,bend);
}

// hexapod lik solver
bool limb_solver1(int limbi, extvec& pos_limb, extvec& joint_angles, bool bend){
  int ysign = (limbi % 2 == 0)? 1 : -1;
  return limb_solver_yxx(limbi,pos_limb,joint_angles,limb_ls,ysign,bend);
}

// spider lik solver
bool limb_solver2(int limbi, extvec& pos_limb, extvec& joint_angles, bool bend){
  int ysign = (limbi % 2 == 0)? 1 : -1;
  return limb_solver_zxx(limbi,pos_limb,joint_angles,limb_ls1,ysign,bend);
}

// Solves for 3-link limb knew bend with hinge axis y-x-x.
bool bend_solver_yxx(int limbi, extvec& pos, const extvec& angles, const double* ls, int ysign, bool pos_flag){
  int s0 = ysign;
  double alpha0, alpha1, alpha2;
  angles.get_components(alpha0,alpha1,alpha2);

  if(pos_flag){
    double l0=ls[0], l1=ls[1], l2=ls[2];
    int s1 = (alpha2 < 0) ? 1 : -1;

    double ca = cos(alpha2);
    double l = sqrt(l1*l1+l2*l2+2*l1*l2*ca);
    double theta = s1*acos((l1+l2*ca)/l) - alpha1 + (1-s0)*M_PI/2;
    double phi = - alpha0;

    // rotated frame, x->x1, y->z1, z->-y1 (?)
    double stl = sin(theta)*l;
    double x1 = stl*sin(phi);
    double y1 = stl*cos(phi);
    double z1 = s0*l0 + cos(theta)*l;

    pos.set(x1,y1,z1);
  }

  return (alpha2*s0 < 0);
}

// quad bend solver
bool bend_solver0(int limbi, extvec& pos, extvec& angles, bool pos_flag){
  int ysign = (limbi<2)? 1 : -1;
  return bend_solver_yxx(limbi,pos,angles,limb_ls,ysign,pos_flag);
}

// hex bend solver
bool bend_solver1(int limbi, extvec& pos, extvec& angles, bool pos_flag){
  int ysign = (limbi % 2 == 0)? 1 : -1;
  return bend_solver_yxx(limbi,pos,angles,limb_ls,ysign,pos_flag);
}

// spider bend solver
// Solver is unfinished! Always returns true bend.
// TODO: finish bend solver.
bool bend_solver2(int limbi, extvec& pos, extvec& angles, bool pos_flag){
  return true;
}


// Child is limb's top-link node. 
liklimb::liklimb(int limbi_, const modelnode* child_){
  limbi = limbi_;
  child = child_;
  parent = child->get_parent();
  setup_joint_values();
  limb_bend = true;
}

// Pulls value refeerences from model nodes to values array.
void liklimb::setup_joint_values(){
  const modelnode* mnode = child;
  for(int i=0;i<3;i++){
    double* val = mnode->get_joint()->get_values();
    values[i] = val;
    mnode = mnode->get_first_child();
  }
}

// Solves for joint values given foot position.
// If successful, sets values in the model.
// Otherwise, program is aborted.
void liklimb::place_limb(const extvec& pos_ground){
  extvec pos_limb;
  poslimb(pos_ground,pos_limb);

  extvec joint_angles;
  if(!solver_func[0](limbi,pos_limb,joint_angles,limb_bend)){
    cout << "limbi = " << limbi << endl;
    cout << "pos_ground: ";
    pos_ground.print();
    extvec pos0;
    get_pos0(pos0);
    cout << "pos0: ";
    pos0.print();
    exit(1);
  }
  
  //joint_angles.print();
  set_joint_values(joint_angles);

  //pos_ground.print();pos_limb.print();exit(1);
}

// Computes foot position in the hip's joint frame.
void liklimb::poslimb(const extvec& pos_ground, extvec& pos_limb){
  modeljoint* joint = child->get_joint();
  joint->compute_A_ground(parent->get_A_ground());
  affine A (*joint->get_A_ground());
  A.invert_rigidbody();
  A.mult(pos_ground,pos_limb);
}

// Sets joint values in the model.
void liklimb::set_joint_values(const extvec& joint_values){
  const double* p = joint_values.get_data();
  double** p1 = values;
  for (int i=0;i<3;i++) {**p1++ = *p++;}
}

// Gets hip position (using child body frame's position).
void liklimb::get_pos0(extvec& pos){
  const affine* A = child->get_A_ground();
  A->get_translation(pos);
}

// Gets foot model node.
modelnode* liklimb::get_foot(){
  return child->get_first_child()->get_first_child(); // limb is assumed to have three links
}

// Tests LIK solver functions.
// See liksolver::solver_test() for more info.
// Test is only suitable for yxx limbs.
void liklimb::solver_test_yxx(int n){
  int m = 1; // use m = 1(2) for testing bend and pos (pos only)
  for(int i=0;i<n;i++){
    extvec angles0, pos, angles1, angles2, angles[2];
    for(int j=0;j<3;j++){
      angles0.set_v(j,(2*randf()-1)*M_PI);
    }
    bool bend = solver_func[1](limbi,pos,angles0,true);
    bool skip_flag = false;
    double score = 1;
    for(int k=0;k<m;k++){
      bool bendk = bool((int(bend)+k)%2);
      solver_func[0](limbi,pos,angles[k],bendk);
      angles[k].subtract(angles0);
      double* p = angles[k].get_data();
      for(int l=0;l<3;l++){mod_twopi(*p++);}
      double d = angles[k].get_v(0)-M_PI;
      mod_twopi(d);
      if(fabs(d)<1e-3){skip_flag = true;} // phi+M_PI branch is ignored
      score *= angles[k].norm();
    }
    if(skip_flag){continue;}

    if(score>1e-3){cout<<"Test fail:"<<endl;
      cout<<"i = "<<i<<endl<<"limbi = "<<limbi<<endl;
      cout<<"angles0: ";angles0.print();
      cout<<"pos: ";pos.print();
      for(int k=0;k<m;k++){
	cout<<"dangles "<<k<<": ";angles[k].print();
      }
      exit(1);
    }
  }
}

bool liklimb::bend_from_angles(extvec& angles){
  extvec pos;
  return solver_func[1](limbi, pos, angles, false);
}


/*void kinematicmodel::set_lik(){
  lik = new liksolver (this);
  }*/

