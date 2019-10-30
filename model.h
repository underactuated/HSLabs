#ifndef MODEL_H
#define MODEL_H

#include <list>
#include <set>

#include "matrix.h"
#include "visualization.h"

enum joint_type {free6, hinge, slider};

class modeljoint{
  joint_type type;
  affine A_parent; // parent-to-joint A
  affine A_ground; // ground-to-joint A
  double* values;
  int val_size;
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
};


class modelnode{
  affine A_tobody; // parent-to-body A if no joint or join-to-body A if joint
  affine A_ground; // ground-to-body A
  modeljoint* joint;
  list<modelnode*> child_nodes;
  modelnode* parent;
public:
  modelnode(const double* pos, const affine* A);
  ~modelnode();
  inline affine* get_A_ground(){return &A_ground;}
  inline affine* get_A_tobody(){return &A_tobody;}
  inline modeljoint* get_joint(){return joint;}
  inline modelnode* get_first_child(){return child_nodes.front();}
  inline modelnode* get_parent(){return parent;}
  void add_child(modelnode* child);
  void print();
  list<double*> make_joint(joint_type type, const double* pos, const double* axis);
  list<double*> make_joint(joint_type type, const double* pos_);
  void recompute_A_ground(const affine& A);
private:
  void compute_A_ground(const affine* A);
};

class liksolver;

class kinematicmodel{
  modelnode* mrootnode;
  visualizer* vis;
  vector<odepart*> odeparts;
  vector<double*> joint_values;
  string xmlfname;
  liksolver* lik;
  bool vis_flag;
  vector<modelnode*> mnodes;
public:
  kinematicmodel(bool vis_flag);
  ~kinematicmodel();
  inline modelnode* get_mnode(int i){return mnodes[i];} // gets model node by index (print model to see indexing)
  inline string get_xmlfname() const {return xmlfname;}
  inline vector<double*>* get_joint_values(){return &joint_values;}
  inline liksolver* get_lik(){return lik;}
  inline visualizer* get_vis(){return vis;}
  inline vector<odepart*>* get_odeparts(){return &odeparts;}
//inline const vector<odepart*>* get_odeparts() const {return &odeparts;}
  inline int number_of_motor_joints() const {return joint_values.size()-6;}
  inline bool if_loaded(){return !(xmlfname == "");}
  inline odepart* get_odepart(int i){return odeparts[i];}
  inline bool if_vis(){return vis_flag;}
  void set_vis();
  void load_fromxml(string fname);
  modelnode* mnode_from_xnode(xml_node<>* xnode, const affine* A_parent);
  void make_odepart(const xml_node<>* xnode, modelnode* mnode);
  void make_joint(const xml_node<>* xnode, modelnode* mnode);
  void orient_bodys();
  void draw();
  void recompute_modelnodes();
  void print();
  void print(int datail_level);
  int get_config_dim() const;
  void set_jvalues_with_lik(const double* rec);
  void set_jvalues(const double* values);
  void get_jvalues(double* values);
  void set_ode_joints();
  void orient_torso(const extvec* orientation);
  void get_foot_mnodes(set<modelnode*>& foot_set);
private:
  void set_lik();
};

#endif
