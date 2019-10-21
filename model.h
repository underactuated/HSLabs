#include <list>
#include <set>

enum joint_type {free6, hinge, slider};

struct modeljoint{
  joint_type type;
  affine A_parent; // parent-to-joint A
  affine A_ground; // ground-to-joint A
  double* values;
  int val_size;
  modeljoint(joint_type type_, affine& A);
  ~modeljoint();
  list<double*> get_values_list();
  double* get_values(){return values;}
  void transformation(affine& A);
  void compute_A_ground(affine* A);
  affine* get_A_parent(){return &A_parent;}
  affine* get_A_ground(){return &A_ground;}
  joint_type get_type(){return type;}
};


struct modelnode{
  affine A_tobody; // parent-to-body A if no joint or join-to-body A if joint
  affine A_ground; // ground-to-body A
  modeljoint* joint;
  list<modelnode*> child_nodes;
  modelnode* parent;
  modelnode(double* pos, affine* A);
  ~modelnode();
  void add_child(modelnode* child);
  affine* get_A_ground(){return &A_ground;}
  affine* get_A_tobody(){return &A_tobody;}
  void print();
  list<double*> make_joint(joint_type type, double* pos, double* axis);
  list<double*> make_joint(joint_type type, double* pos_);
  void recompute_A_ground(affine& A);
  modeljoint* get_joint(){return joint;}
  modelnode* get_first_child(){return child_nodes.front();}
  modelnode* get_parent(){return parent;}
private:
  void compute_A_ground(affine* A);
};

struct liksolver;

struct kinematicmodel{
  modelnode* mrootnode;
  visualizer* vis;
  vector<odepart*> odeparts;
  vector<double*> joint_values;
  string xmlfname;
  liksolver* lik;
  bool vis_flag;
  vector<modelnode*> mnodes;
  kinematicmodel(bool vis_flag);
  ~kinematicmodel();
  void set_vis();
  void load_fromxml(string fname);
  modelnode* mnode_from_xnode(xml_node<>* xnode, affine* A_parent);
  void make_odepart(xml_node<>* xnode, modelnode* mnode);
  void make_joint(xml_node<>* xnode, modelnode* mnode);
  void orient_bodys();
  void draw();
  void recompute_modelnodes();
  void print();
  void print(int datail_level);
  modelnode* get_mnode(int i);
  string get_xmlfname(){return xmlfname;}
  vector<double*>* get_joint_values(){return &joint_values;}
  liksolver* get_lik(){return lik;}
  visualizer* get_vis(){return vis;}
  vector<odepart*>* get_odeparts(){return &odeparts;}
  int get_config_dim();
  void set_jvalues_with_lik(double* rec);
  void set_jvalues(double* values);
  void get_jvalues(double* values);
  int number_of_motor_joints(){return joint_values.size()-6;}
  void set_ode_joints();
  void orient_torso(extvec* orientation);
  void get_foot_mnodes(set<modelnode*>& foot_set);
  bool if_loaded(){return !(xmlfname == "");}
private:
  void set_lik();
};

