/////////////////////////////////////////////////
// dynrec.c: declares classes that determine dynamic
// properties of model parts and store various 
// quantitites reflecting part's dynamics along
// the trajectory. These classes are also used
// to encode relationship between forces, torques
// and accelerations of model's parts. We use
// sparse matrix representaion for computational
// efficiency.
/////////////////////////////////////////////////
#ifndef DYNREC_H
#define DYNREC_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "matrix.h"
#include "model.h"
#include "pergen.h"
#include "periodic.h"

using namespace Eigen;

typedef SparseMatrix<double> SpMat; // column-major sparse matrix

// Class dynpart ( = dynamics part) contains geometric, inertial
// and contact properties of model's part. Includes access to
// corresponding ode-part and model-node objects.
class dynpart{
  int id, parent_id; // part and its parent's ids
  extvec com_pos, joint_pos, foot_pos;
  const odepart* opart; // corresponding ode-part
  modelnode* mnode; // corresponding model node
  double mass;
  affine inertia; // TODO: maybe replace affine with dMatrix3?
  bool foot_flag; // part is a foot if foot_flag
public:
  dynpart(const odepart* opart_){opart = opart_;}
  inline extvec* get_com_pos(){return &com_pos;}
  inline extvec* get_joint_pos(){return &joint_pos;}
  inline extvec* get_foot_pos(){return &foot_pos;}
  inline int get_parent_id(){return parent_id;}
  inline double get_mass(){return mass;}
  inline bool if_foot(){return foot_flag;}
  void setup(map<modelnode*,int>& mnode_id_map);
  void print();
  void recompute();
  const affine* get_A_ground();
  affine* get_inertia_tensor();
  void setup_foot(const set<modelnode*>& foot_set);
  void get_joint_zaxis(extvec& axis);
private:
  void set_joint_pos();
  void set_com_pos();
  void set_inertial_params();
  void set_foot_pos();
};

// Class dynrecord ( = dynamics record) computes and stores
// various quantities of the model parts pertinent to dynamic
// properties of the trajectory (for a single time moment).  
// It also constructs a linear system of equations describing
// relationship between the joint forces and torques and 
// the model parts forces and torques.
// We denote forces and torques applied at joints (including
// free and fixed joints) and contacts x and call it joint-ft
// vector; We denote forces and torques applied at part COMs
// f and call it com-ft vector; we use B to denoted the matrix
// that connects x and f, and we call if ft matrix ( = force-
// torque matrix): B * x = f. If the model has n parts (each
// part contains exactly one joint) and k contacts, then 
// dimensions of x, f and B are: dim(x) = 6n + 3k, dim(f) = 6n
// and dim(B) = 6n x (6n+3k). Note that a motor torque is
// the projection of the full joint torque on the joint axis.
// TODO: maybe store data in a more compact/efficient represention
class dynrecord{
  int n; // number of parts
  int nf; // mumber of feet
  extvec *pos, *jpos, *vel, *mom, *mom_rate, *acc, *ust, *ang_vel, *ang_mom, *ang_mom_rate, *fpos, *jzaxis; // COM quantities: pos = position, vel = velocity, mom = momentum, mom_rate = time derivative of mom, acc = acceleration, ust = u*sin(theta) where u is the unit axis vectora, ang_vel = angular velocity, ang_mom = angular momentum, ang_mom_rate = time derivative of ang_mom; joint quantities: jpos = position, jzaxis = z component of joint axis; fpos = foot position.
  affine *rot; // part rotations
  bool* contacts; // foot contact states
public:
  dynrecord(int n, int nf);
  ~dynrecord();
  inline extvec* get_mom_rate(){return mom_rate;}
  inline extvec* get_ang_mom_rate(){return ang_mom_rate;}
  inline extvec* get_pos(){return pos;}
  inline extvec* get_jpos(){return jpos;}
  inline extvec* get_jzaxis(){return jzaxis;}
  void initialize(vector<dynpart*>& dynparts, double rcap);
  void print(int id);
  void compute_ders(int stage, const dynrecord* prev_rec, const dynrecord* next_rec, double dt, periodic* per);
  void set_forcetorque_system(SpMat& B, VectorXd& f, const int* parentis, const double* masses);
  int get_ncontacts();
  void set_forcetorque_system_contacts(SpMat& B, const int* footis);
  void set_forcetorque_system_contacts(SpMat& B, const int* footis, bool contact_feet_flag);
private:
  void compute_ders(extvec *p_der, const extvec *p_func_prev, const extvec *p_func_next, double dt);
  void compute_ang_mom(periodic* per);
  void compute_mom(const double* masses);
  void ftsys_forces(int i, int pi, SpMat& B, VectorXd& f);
  void ftsys_torques(int i, int pi, SpMat& B, VectorXd& f);
  void ftsys_gravity(int i, VectorXd& f, const double* ms);
  void ftsys_contact_forces(int i, int ci, SpMat& B);
  void ftsys_contact_torques(int fi, int i, int ci, SpMat& B);
};

#endif
