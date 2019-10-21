///////////////////////////////////////////////////////////
// lik design outline from todo list, to be improved
// --------------------------------------------------------
// lik outline:
// LIK = Limb Inverse Kinematics
// lik has array of liklimbs objects (for each limb). Each liklimb contains first
// leg-link model node and its parent. Lik function takes as input: limb index and
// foot tip position; it sets joint values to the limb; optionally, it can use model
// node to verify correctness of its output.
///////////////////////////////////////////////////////////

typedef bool(*SolverFuncType)(int, extvec&, extvec&, bool);

struct liklimb;

// for now we assume that limbs are identical
struct liksolver{
  int index;
  vector<liklimb*> limbs;
  liksolver(kinematicmodel* model);
  SolverFuncType* solver_func;
  ~liksolver();
  void place_limb(int limbi, double x, double y, double z);
  void place_limbs(double* rec);
  void get_limb_pos0(int limbi, extvec& pos);
  void print_limb_pos0s();
  int get_number_of_limbs(){return limbs.size();}
  vector<liklimb*>* get_limbs(){return &limbs;}
  void solver_test(int n);
private:
  void set_limbs(kinematicmodel* model);
};

struct liklimb{
  int limbi;
  modelnode *parent, *child;
  double* values[3];
  SolverFuncType* solver_func;
  bool limb_bend;
  liklimb(int limbi_, modelnode* child_);
  void setup_joint_values();
  void place_limb(extvec& pos_ground);
  void set_joint_values(extvec& joint_values);
  void get_pos0(extvec& pos);
  void set_solver_func(SolverFuncType* f){solver_func = f;}
  modelnode* get_foot();
  void solver_test_yxx(int n);
  bool bend_from_angles(extvec& angles);
  void set_bend(bool bend){limb_bend = bend;}
private:
  void poslimb(extvec& pos_ground, extvec& pos_limb);
};


