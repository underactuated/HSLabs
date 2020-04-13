#include "model.h"
#include "lik.h"

// A_parent (arg A), of joint host node's parent,
// is supplied at construction
modeljoint::modeljoint(joint_type type_, const affine& A){
  type = type_;
  //A_parent.copy(A);
  A_parent = A;
  switch(type){
  case free6: val_size = 6; break;
  default: val_size = 1; break;
  }
  values = new double [val_size];
  double* p = values;
  for(int i=0;i<val_size;i++){*p++ = 0;}
}

modeljoint::~modeljoint(){
  delete [] values;
}

// Returns joint values as a list.
list<double*> modeljoint::get_values_list(){
  double* p = values;
  list<double*> values_list;
  for(int i=0;i<val_size;i++){
    values_list.push_back(p);
    p++;
  }
  return values_list;
}

// joint transformation corresponding to its value
// free6: rotation around origin, followed by translation
// TODO: implement slider joint
void modeljoint::transformation(affine& A){
  switch(type){
  case free6: {
    extvec pos;
    dMatrix3 rot;
    double* p = values;
    pos.set(p);
    p = values+3;
    dRFromEulerAngles(rot,*p,*(p+1),*(p+2));
    A.set_rotation(rot);
    A.translate(pos);
  } break; 
  case hinge: {
    double val = *values;
    double c = cos(val), s = sin(val);
    A.set_unity();
    for(int i=0;i<2;i++){
      A.set_a(i,i,c);
      A.set_a(i,1-i,(2*i-1)*s);
    }
  } break;
  default:
    cout<<"ERROR: joint type "<<type<<" transformation has not been implemented"<<endl;exit(1);
    break;
  }
}

// Computes joint's A_ground from parent’s A_ground (arg A).
void modeljoint::compute_A_ground(const affine* A){
  A->mult(A_parent,A_ground);
}

void modeljoint::print(){
  cout << "--- model joint ---" << endl;
  cout << "type = " << type << endl;
  cout << "A_parent:" << endl;
  A_parent.print();
  cout << "A_ground:" << endl;
  A_ground.print();
}


// Position pos, relative to parent, and parent's 
// A_ground (arg A) are supplied at construction
modelnode::modelnode(const double* pos, const affine* A){
  A_pj_body.set_translation(pos);
  compute_A_ground(A);
  joint = NULL;
  parent = NULL;
}

modelnode::~modelnode(){
  if (joint) {delete joint;}
}

void modelnode::add_child(modelnode* child){
  child_nodes.push_back(child);
  child->parent = this;
}

// computes A_ground at node construction only (before a joint is added)
// from parent's A_ground (arg A)
void modelnode::compute_A_ground(const affine* A){
  //A_ground.copy(*A);
  A_ground = (*A);
  A_ground.mult(A_pj_body);
}

void modelnode::print(int detail_level) const {
  cout << "--- model node ---" << endl;
  cout << "A_pj_body:" << endl;
  A_pj_body.print();
  cout << "A_ground:" << endl;
  A_ground.print();
  cout << "has " << child_nodes.size() << " child nodes" << endl;
  if(detail_level > 0){
    if(joint != NULL){joint->print();}
  }
}

// Makes a joint specified by an axis, such as hinge or slider,
// (the latter remains to be implemented).
list<double*> modelnode::make_joint(joint_type type, const double* pos_, const double* axis){
  switch (type) {
  case hinge: {
    extvec pos, v;
    pos.set(pos_);
    v.set(axis);
    
    dMatrix3 rot;
    rot_ztov(rot,v);

    affine A1;
    // setting joint
    A1.set_rotation(rot);
    affine A (A_pj_body);
    A.mult(A1);
    A.translate(pos);
    joint = new modeljoint (type,A);

    // modifying A_pj_body
    dVector3 pos_dv;
    pos.to_dvec(pos_dv);
    affine_from_posrot(A_pj_body,pos_dv,rot);
    A_pj_body.invert_rigidbody();
  } break;
  default:
    cout<<"ERROR: "<<type<<" not implemented"<<endl;exit(1);
    break;
  }
  return joint->get_values_list();
}

// Makes a free joint.
list<double*> modelnode::make_joint(joint_type type, const double* pos_){
  switch (type) {
  case free6: {
    extvec pos;
    pos.set(pos_);

    // setting joint
    affine A (A_pj_body);
    A.translate(pos);
    joint = new modeljoint (type,A);
    
    // modifying A_pj_body
    A.set_unity();
    pos.times(-1);
    A.translate(pos);
    //A_pj_body.copy(A);
    A_pj_body = A;
  } break;
  default:
    cout<<"ERROR: "<<type<<" not implemented"<<endl;exit(1);
    break;
  }
  return joint->get_values_list();
}

// relation between A_ground of two bodies connected by a joint:
// let mnode1 be parent of mnode2 connected by joint
// let B = mnode1.A_ground, C = mnode2.A_ground, D = joint.A_parent
// E = joint.transformation, F = mnode2.A_pj_body
// then C = B * D * E * F
//////////////////////////
// Recomputes A_ground from parent's A_ground
void modelnode::recompute_A_ground(const affine& A){
  if(joint){
    joint->compute_A_ground(&A);
    //A_ground.copy(*joint->get_A_ground());
    A_ground = (*joint->get_A_ground());
    affine A_joint_pastjoint;
    joint->transformation(A_joint_pastjoint);
    A_ground.mult(A_joint_pastjoint);
    A_ground.mult(A_pj_body);
  } else {
    //A_ground.copy(A);
    A_ground = A;
    A_ground.mult(A_pj_body);
  }
  list<modelnode*>::iterator it = child_nodes.begin();
  for(;it!=child_nodes.end();it++){
    (*it)->recompute_A_ground(A_ground);
  }
}


// If !vis_flag, the model is not represented in a visualizer.
// This is used in ghost walking, as a ghost model does
// not have a physical presence in ODE world.
kinematicmodel::kinematicmodel(bool vis_flag_){
  vis_flag = vis_flag_;
  if(vis_flag){set_vis();}
  mrootnode = NULL;
  lik = NULL;
}

kinematicmodel::~kinematicmodel(){
  if (mrootnode) {delete mrootnode;}
  if (lik) {delete lik;}
  while(odeparts.size()){
    delete odeparts.back();
    odeparts.pop_back();
  }
}

// Loads a model from xml file.
void kinematicmodel::load_fromxml(string fname){
  xmlfname = fname;
  file<> xmlFile(fname.c_str());
  xml_document<> doc;
  doc.parse<0>(xmlFile.data());

  xml_node<> *node = doc.first_node("mujoco");
  if(!node){cout << "ERROR: not a mujoco file" << endl; exit(1);}
  node = node->first_node("worldbody");
  node = node->first_node("body");
  affine A;
  A.set_unity();
  mrootnode = mnode_from_xnode(node,&A);
  lik = new liksolver (this);
  recompute_modelnodes();
  // Recomputation is needed for correct joint.A_ground,
  // as can be checked by printing model->print(5)
  // before and after recomputation.
}

// Makes model node according to xml node specifications.
// Arg A is parent's A_ground
modelnode* kinematicmodel::mnode_from_xnode(xml_node<>* xnode, const affine* A){
  double pos[3];
  xmlnode_attr_to_val(xnode,"pos",pos);
  modelnode* mnode = new modelnode (pos,A);
  make_odepart(xnode,mnode);
  make_joint(xnode,mnode);
  xml_node<>* xchild = xnode->first_node("body");
  while(xchild){
    //cout << xchild->name() << endl;
    modelnode* mchild = mnode_from_xnode(xchild,mnode->get_A_ground());
    xchild = xchild->next_sibling("body");
    mnode->add_child(mchild);
  }
  return mnode;
}

// Makes odepart if vis_flag.
void kinematicmodel::make_odepart(const xml_node<>* xnode, modelnode* mnode){
  mnodes.push_back(mnode);
  if(!vis_flag){return;}
  odepart* part = new odepart;
  part->make(xnode,mnode,vis);
  odeparts.push_back(part);
}

// Makes model joint for model node mnode.
void kinematicmodel::make_joint(const xml_node<>* xnode, modelnode* mnode){
  xml_node<>* joint_node = xnode->first_node("joint");
  if(!joint_node){return;}
  string type = joint_node->first_attribute("type")->value();
  //cout << "joint type: " << type << endl;
  list<double*> values;
  if(type == "free"){
    double pos[3];
    xmlnode_attr_to_val(joint_node,"pos",pos);
    values = mnode->make_joint(free6,pos);
  } else if(type == "hinge"){
    double pos[3], axis[3];
    xmlnode_attr_to_val(joint_node,"pos",pos);
    xmlnode_attr_to_val(joint_node,"axis",axis);
    values = mnode->make_joint(hinge,pos,axis);
  }
  joint_values.insert(joint_values.end(),values.begin(),values.end());
}

// Orients ode bodys according to model nodes.
void kinematicmodel::orient_odebodys(){
  vector<odepart*>::iterator it = odeparts.begin();
  for(;it!=odeparts.end();it++){
    odepart* part = (*it);
    dVector3 pos;
    dMatrix3 rot;
    //part->get_body_posrot_from_frame(pos,rot);
    part->get_odebody_posrot_from_body(pos,rot);
    dBodyID odebody = part->get_odebody();
    dReal* p = pos;
    dBodySetPosition(odebody,*p,*(p+1),*(p+2));
    dBodySetRotation(odebody,rot);
  }
}

// Launches visualizer's simulation/drawing loop.
void kinematicmodel::draw(){
  vis->draw();
}

// Recomputes kinematic model tree. Should be used
// after configuration values are changed.
void kinematicmodel::recompute_modelnodes() const {
  affine A;
  A.set_unity();
  mrootnode->recompute_A_ground(A);
}

void kinematicmodel::print() const {
  print(0);
}

void kinematicmodel::print(int detail_level) const {
  cout << "--- kinematic model ---" << endl;
  cout << "loaded from " << xmlfname << endl;
  int config_dim = get_config_dim();
  cout << "number of DoFs = " << config_dim <<endl;
  cout << "configuration:" << endl;
  for(int i=0;i<config_dim;i++){
    if(i){cout << " ";}
    cout << *joint_values[i];
  }
  cout << endl;
  int odeparts_size = odeparts.size();
  cout << "number of ode parts: " << odeparts_size << endl;
  cout << "ode parts:" << endl;
  for(int i=0;i<odeparts_size;i++){
    if(detail_level == 0){
      cout << i << ": " << odeparts[i]->get_part_name() << endl;
    }
    if(detail_level > 0){odeparts[i]->print(detail_level-1);}
  }
}

// Returns configuration dimensionality,
// assumed to be the total number of joint values.
int kinematicmodel::get_config_dim() const {
  return joint_values.size();
}

// Sets joint values with LIK, via torso orientation
// and limb positions, contained in rec ( = record).
void kinematicmodel::set_jvalues_with_lik(const double* rec) const {
  const double* p = rec;
  for(int i=0;i<6;i++){*joint_values[i] = *p++;}
  recompute_modelnodes();
  lik->place_limbs(p);
}

// Sets joint values directly.
void kinematicmodel::set_jvalues(const double* values) const {
  const double* p = values;
  vector<double*>::const_iterator it = joint_values.begin();
  for(;it!=joint_values.end();it++){*(*it) = *p++;}
}

void kinematicmodel::get_jvalues(double* values) const {
  double* p = values;
  vector<double*>::const_iterator it = joint_values.begin();
  for(;it!=joint_values.end();it++){*p++ = *(*it);}
}

// ODE joints are set after kinematic model is completed.
void kinematicmodel::set_ode_joints(){
  map<modelnode*,odepart*> mnode_part_map;
  vector<odepart*>::iterator it = odeparts.begin();
  for(;it!=odeparts.end();it++){
    odepart* part = (*it);
    mnode_part_map[part->get_mnode()] = part;
  }
  for(it=odeparts.begin();it!=odeparts.end();it++){
    odepart* part = (*it);
    modelnode* mnode = part->get_mnode();
    modelnode* pmnode = mnode->get_parent();
    if(pmnode){
      odepart* ppart = mnode_part_map[pmnode];
      modeljoint* joint = mnode->get_joint();
      if(!joint){part->make_fixed_joint(ppart,vis);}
      else {
	switch (joint->get_type()) {
	case hinge:
	  part->make_hinge_joint(ppart,vis);
	  break;
	default: break;
	}
      }
    } 
  }
}

// Sets new torso irientation and recomputes model.
void kinematicmodel::orient_torso(const extvec* orientation) const {
  for(int i=0;i<2;i++){
    const double* p = orientation[i].get_data();
    for(int j=0;j<3;j++){*joint_values[j+i*3] = *p++;}
  }
  recompute_modelnodes();
}

// Gets foot model nodes.
void kinematicmodel::get_foot_mnodes(set<modelnode*>& foot_set) const {
  const vector<liklimb*>* limbs = get_lik()->get_limbs();
  vector<liklimb*>::const_iterator it = limbs->begin();
  for(;it!=limbs->end();it++){
    foot_set.insert((*it)->get_foot());
  }
}

