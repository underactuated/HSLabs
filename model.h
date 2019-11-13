/////////////////////////////////////////////////
// model.h: declaration of classes related to
// kinematic model implementation. Kinematic model
// can be viewed as a computation engine for
// relating model's configuration coordinated
// (such as torso position, orientation and all
// the joint angles) to cortesian coordinates of
// its parts. It also laverages inverse kinematics
// for converting torso and foot coordinates into
// a model configuration. Kinematic model is fully
// specified in a xml file following URDF-ish format.
// Inverse kinematics requires imlementation of one
// or two functions (limb_solver# and bend_solver#)
// in lik.cpp. The second function is needed for
// extreemly rough terrain stability. 
/////////////////////////////////////////////////
#ifndef MODEL_H
#define MODEL_H

#include <list>
#include <set>

#include "matrix.h"
#include "visualization.h"

enum joint_type {free6, hinge, slider};

// Class modeljoint implements transformation
// associated with joint value(s)
// General rules for affine transformation notations:
// A_frame0_frame1 denotes transformation from
// frame0 to frame1. If frame1 is this object,
// it can be omitted.
class modeljoint{
  joint_type type;
  affine A_parent; // parent-to-joint (transformation) A
  affine A_ground; // ground-to-joint (transformation) A
  double* values;
  int val_size; // number of configuration parameters (values) at the joint
public:
  modeljoint(joint_type type_, const affine& A);
  ~modeljoint();
  inline double* get_values(){return values;}
  inline affine* get_A_parent(){return &A_parent;}
  inline affine* get_A_ground(){return &A_ground;}
  inline joint_type get_type(){return type;}
  list<double*> get_values_list();
  void transformation(affine& A);
  void compute_A_ground(const affine* A);
  void print();
};


// Class modelnode implements kinematic model node
// ( = model node). A model node corresponds to a rigid part
// of a robot model, optionally containing a joint, connecting
// the model node to its parent. No joint in the model node
// (no joint for a body in xml file) corresponds to fixed ODE joint.
// We distinguish between body and ODE body ( = ode-body).
// Body and body frame referr to body in xml file.
// ode-body and ode-body frame referr to ODE body and corresponding
// geom. The two are different in general, though sometimes 
// their frames coniside (when odepart.A_body_geom is an identity).
class modelnode{
  affine A_pj_body; // pj = parent or joint: parent if no joint, joint otherwise
  affine A_ground; // ground-to-body A
  modeljoint* joint;
  list<modelnode*> child_nodes;
  modelnode* parent;
public:
  modelnode(const double* pos, const affine* A);
  ~modelnode();
  inline const affine* get_A_ground() const {return &A_ground;}
  inline affine* get_A_pj_body(){return &A_pj_body;}
  inline modeljoint* get_joint() const {return joint;}
  inline modelnode* get_first_child() const {return child_nodes.front();}
  inline modelnode* get_parent() const {return parent;}
  void add_child(modelnode* child);
  void print(int detail_level);
  void print(){print(0);}
  list<double*> make_joint(joint_type type, const double* pos, const double* axis);
  list<double*> make_joint(joint_type type, const double* pos_);
  void recompute_A_ground(const affine& A);
private:
  void compute_A_ground(const affine* A);
};

class liksolver;

// Class kinematicmodel implements and manipulates kinematic model.
// Kinematic model is represented by a (kinematic) tree,
// with model nodes as tree nodes.
// (Kinematic) model has access to LIK (limb inverse kinematics) solver.
// It has optinal representation in visualizer, in which case
// it stores odeparts. The model can be loaded from xml file.
class kinematicmodel{
  modelnode* mrootnode; // currently, the root node is always torso node 
  visualizer* vis;
  vector<odepart*> odeparts;
  vector<double*> joint_values;
  string xmlfname; // xml file name
  const liksolver* lik; // LIK solver
  bool vis_flag; // visualizer it optional, not used for ghost model
  vector<modelnode*> mnodes;
public:
  kinematicmodel(bool vis_flag);
  ~kinematicmodel();
  inline const modelnode* get_mnode(int i) const {return mnodes[i];} // Gets model node by index (print model to see indexing).
  inline string get_xmlfname() const {return xmlfname;}
  inline vector<double*>* get_joint_values(){return &joint_values;}
  inline const liksolver* get_lik() const {return lik;}
  inline visualizer* get_vis() const {return vis;}
  inline const vector<odepart*>* get_odeparts() const {return &odeparts;}
  inline int number_of_motor_joints() const {return joint_values.size()-6;}
  inline bool if_loaded(){return !(xmlfname == "");}
  inline const odepart* get_odepart(int i) const {return odeparts[i];}
  inline bool if_vis() const {return vis_flag;}
  //void set_vis();
  void load_fromxml(string fname);
  void orient_odebodys();
  void draw();
  void recompute_modelnodes() const;
  void print();
  void print(int detail_level);
  int get_config_dim() const;
  void set_jvalues_with_lik(const double* rec) const;
  void set_jvalues(const double* values) const;
  void get_jvalues(double* values) const;
  void set_ode_joints();
  void orient_torso(const extvec* orientation) const;
  void get_foot_mnodes(set<modelnode*>& foot_set) const;
private:
  void set_vis();
  modelnode* mnode_from_xnode(xml_node<>* xnode, const affine* A);
  void make_odepart(const xml_node<>* xnode, modelnode* mnode);
  void make_joint(const xml_node<>* xnode, modelnode* mnode);
};

#endif
