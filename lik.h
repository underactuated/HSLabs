/////////////////////////////////////////////////
// lik.h: declaration of classes and functions
// needed for inverse kinematics computation (IK).
// Their primary purpose is to compute model's 
// configuraiton (joint angles), that realize 
// a desired arrangement of foot poisitions and
// torso orientation.
// Currently, the limbs are interchangible, so 
// they have a shared set of IK functions.
// IK is primarily used for generating walking cycle
// trajectories and in ghost walking transformation.
/////////////////////////////////////////////////
#ifndef LIK_H
#define LIK_H

#include "matrix.h"
#include "model.h"

typedef bool(*SolverFuncType)(int, extvec&, extvec&, bool);

class liklimb;

// Class liksolver (for Limb inverse kinematics solver,
// or LIK solver) solves IK by setting up and calling
// liklimb objects, responsible for solving IK for
// every limb. limb_solver function -- the central 
// element in computing limb IK -- defined in the source
// file, solves IK in hip's joint frame for a given knee
// bend. (There are two branches of solution, realizing
// a desired foot position, with positive and negative
// knee bends.) Additionally, we define bend_solver
// function that determines bend given joint angles.
// It is needed for ghost walking, because the ghost
// is enforced to have the same bends as the model.
// For now we assume identical limbs, but in principle
// limb differentiation is supported by defining
// liklimb object separately for each limb.
// liksolver also implements (so far partially)
// testing of limb_solver and bend_solver.
class liksolver{
  int index; // enumerates variants for lik solvers. If solver undefined, index = -1.
  vector<liklimb*> limbs;
  SolverFuncType* solver_func;
  double rcap; // foot capsule cap radius
public:
  liksolver(const kinematicmodel* model);
  ~liksolver();
  inline int get_number_of_limbs() const {return limbs.size();}
  inline const vector<liklimb*>* get_limbs() const {return &limbs;}
  inline double get_rcap() const {return rcap;}
  void place_limb(int limbi, double x, double y, double z) const;
  void place_limbs(const double* rec) const;
  void get_limb_hip_pos(int limbi, extvec& pos) const;
  void print_hip_poss() const;
  void solver_test(int n) const;
private:
  void set_limbs(const kinematicmodel* model);
  void set_rcap(const kinematicmodel* model, const vector<int>& limb_inds);
};

///////////////////////////////////////////////////////////
// Class liklimb stores information necessary for solving IK
// for a given limb, indexed by limbi. It includes the first
// leg-link model node (thigh) and its parent, model node joint
// values and limb_ and bend_ solver functions.
//////////////////////////////////////////////////////
class liklimb{
  int limbi;
  const modelnode *parent, *child;
  double* values[3];
  SolverFuncType* solver_func;
  bool limb_bend;
public:
  liklimb(int limbi_, const modelnode* child_);
  inline void set_bend(bool bend){limb_bend = bend;}
  inline void set_solver_func(SolverFuncType* f){solver_func = f;}
  void place_limb(const extvec& pos_ground);
  void set_joint_values(const extvec& joint_values);
  void get_hip_pos(extvec& pos);
  modelnode* get_foot();
  void solver_test_yxx(int n);
  bool bend_from_angles(extvec& angles);
private:
  void setup_joint_values();
  void poslimb(const extvec& pos_ground, extvec& pos_limb);
};

#endif
