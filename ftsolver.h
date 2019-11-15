/////////////////////////////////////////////////
// ftsolver.h: declaration of class forcetorquesolver
// that solves linear force-torque system in various
// settings. In general, the system is underdetermined,
// which requires intruduction of additional constraints.
// We use a perturbation theory approach to ensure,
// that whenver possible, the solver finds a physically
// realizable solution for a given trajectory, as should
// always be the case for a hexapod sufficiently close
// to a (stable) static limit.
/////////////////////////////////////////////////
#ifndef FTSOLVER_H
#define FTSOLVER_H

#include <list>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "periodic.h"
#include "dynrec.h"

using namespace Eigen;

typedef SparseQR<SpMat, COLAMDOrdering<int> > SpSolver; // sparse solver

class dynrecord;

// Class forcetorquesolver for unknown forces and torques.
// A (possibly) underactuated problem is turned into
// a fully actuated problem by allowing non-zero forces
// and torques on the torso (i.e. free joint). Then one
// attempts to find a physically realizable solution
// using perturbation theory by imposing torso penalty
// on forces and torques in the 0th order, while imposing
// motor joint torque penalties in the 1st order. 
class forcetorquesolver{
  int n; // number of dynparts
  int nf; // number of feet
  int *parentis, *footis; // see periodic.h
  double* masses;
  list<pair<int,int> > penal_mask0, penal_mask1; // 0th and 1st order masks
  int mask_l0, mask_l1; // mask lengths
  dynrecord* dynrec; // curent dynamics record
  double* jzaxis; // joint z axis, all dynparts
  VectorXd fts; // last force-torques solution
  list<int> jpart_ids;
public:
  forcetorquesolver(const periodic* per);
  ~forcetorquesolver();
  inline double* get_jzaxis(){return jzaxis;}
  inline VectorXd* get_fts(){return &fts;}
  void print();
  void solve_forcetorques(dynrecord* dynrec_);
  void solve_forcetorques(dynrecord* dynrec_, VectorXd& x, VectorXd& y);
  void switch_torso_penalty(bool force_flag, bool torque_flag);
  void set_action_penalties(VectorXd& c);
  void solve_forces(dynrecord* dynrec_, VectorXd& z, VectorXd& y);
private:
  void solve_forcetorques_particular(VectorXd& x, SpMat& B, VectorXd& f);
  void solve_forcetorques_null_space(SpMat& B, MatrixXd& N);
  void solve_contact_forces(VectorXd& x, VectorXd& y, const MatrixXd& N);
  void extract_N_contact(const MatrixXd& N, MatrixXd& N_cont);
  void set_penal_mask1();
  void set_jzaxis();
  void set_B_and_f_for_zero_cfs(SpMat& B, VectorXd& f);
  void add_torque_constraints_to_B(SpMat& B, VectorXd& f, const VectorXd& z);
  void set_dynrec(dynrecord* dynrec_);
};

#endif
