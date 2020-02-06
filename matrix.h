/////////////////////////////////////////////////
// matrix.h: declaration of augmented matrices and 
// vectors (affine and extvec class correspondingly)
// needed for computing basic chain kinematics.
/////////////////////////////////////////////////
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <ode/ode.h>

using namespace std;

class extvec;

// Class affine represents a 4 x 4 matrix A = [rot,transl;03^T,1]
// cooresponding to a coordinate transformation between two frames
// F and F', where F' can be viewed as the result of rotation
// followed by translation of F by rot and transl respectively.
// (In more details: rot is 3 x 3 rotation matrix, transl is 3-dim 
// translation vector, ^T means transposition, and 03 = [0;0;0].)
// Such transformation is called a rigid body transformation. 
// Then (extended) coordinates r and r' (in F and F' rescpectively)
// of some point are related via r = A * r'. Composition of
// several rigid body transformations A1, A2 ..., is computed
// via matrix multiplication A1 * A2 * ...
class affine{ 
  double a[16]; // matrix is stored column-wise
public:
  inline double* get_data(){return a;}
  inline const double* get_data() const {return a;}
  void set_zeros();
  void set_unity();
  void print_all() const;
  void print() const;
  void set_translation(double x, double y, double z);
  void set_translation(const double* v);
  void copy(const affine& b);
  void copy_transposed(const affine& b);
  void mult(const affine& b);
  void mult(const affine& b, affine& c) const;
  void set_rotation(const dMatrix3& rot);
  void set_rotation(const affine& rot);
  void translate(const extvec& t);
  void transpose();
  void set_a(int i, int j, double val);
  double get_a(int i, int j) const;
  void mult(const extvec& v, extvec& u) const;
  void subtract(const affine& b);
  double norm() const;
  void invert_rigidbody();
  void get_translation(extvec& t) const;
private:
  void print_rows(int n) const;
};

// Class extvec represents a 4-dimensional vector [r;1] 
// corresponding to a 3-dimensional vector r. It is augmented 
// (or extended, hence extvec = extended vector) with 1 as 
// the fourth component.
class extvec{
  double v[4];
public:
  extvec(){v[3] = 1;}
  extvec(double x, double y, double z);
  inline void set_v(int i, double val){v[i] = val;}
  inline double* get_data() {return v;}
  inline const double* get_data() const {return v;}
  inline double get_v(int i) const {return v[i];}
  void set_zeros();
  void print() const;
  void set(double x, double y, double z);
  void set(const double* a);
  void copy(const extvec& u);
  void copy3(const extvec& u);
  void subtract(const extvec& u);
  double norm() const;
  void get_components (double& x, double& y, double& z) const;
  void get_components (double* p) const;
  void cross(const extvec& u, extvec& w) const;
  void cross(const extvec& u);
  void times(double f);
  void normalize();
  double dot(const extvec& u) const;
  void to_dvec(dVector3& u) const;
  void add(const extvec& u);
  double distance(const extvec& u) const;
};

#endif
