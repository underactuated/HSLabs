/////////////////////////////////////////////////
// likredund.h: declaration of redundof class and
// various functions needed for implementing 
// inverse kinematics computation for limbs with
// redundant degrees of freedom (more than 3).
// The primary purpose was to extend the original
// 3DoF design to accomodate Weaver hexapod.
/////////////////////////////////////////////////
#ifndef LIKREDUND_H
#define LIKREDUND_H

void set_ignore_reach_flag_redund(bool value);

bool limb_solver_yxxx(const double* constrs, double* jangles, const double* ls, int ysign, bool bend);
bool limb_solver_zyxxx(const double* constrs, double* jangles, const double* ls, int ysign, bool bend);
bool limb_solver_zyxxx1(const double* constrs, double* jangles, const double* ls, int ysign, bool bend);

typedef void(*RedundFuncType)(int, double*);
typedef void(*LimbRedundFuncType)(affine&, const double*, double*);
typedef void(*RedundGostFuncType)(int, odepart* opart, double*);

class redundof{
  int n, dof;
  RedundFuncType redund_func;
  LimbRedundFuncType limbredund_func;
  RedundGostFuncType redund_gost_func;
public:
  redundof(int index);
  inline int get_dof() const {return dof;}
  inline void set_limb_number(int n_){n = n_;}
  void set_rec_redundof(double* rec) const;
  void constrlimb(affine& A, const double* rec, double* constr);
  void set_rec_redundof_gost(double* rec, vector<odepart*>& foot_oparts) const;
};

void vec_to_direction(const extvec& vec, double* direct);
void direction_to_vec(const double* direct, extvec& vec);
void set_xyz_rotation(affine& A, int axisi, double angle);
void chain_rotation(affine& A, int* axisi, double* q, int n);

#endif
