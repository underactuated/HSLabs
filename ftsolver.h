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

class forcetorquesolver{
  int n; // number of dynparts
  int nf; // number of feet
  int *parentis, *footis;
  double* masses;
  list<pair<int,int> > penal_mask0, penal_mask1;
  int mask_l0, mask_l1;
  dynrecord* dynrec;
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
