#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <ode/ode.h>

using namespace std;

class extvec;

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
  double get_a(int i, int j);
  void mult(const extvec& v, extvec& u) const;
  void subtract(const affine& b);
  double norm() const;
  void invert_rigidbody();
  void get_translation(extvec& t) const;
private:
  void print_rows(int n) const;
};


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
};

#endif
