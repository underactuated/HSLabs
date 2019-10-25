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
  void set_zeros();
  void set_unity();
  void print_all();
  void print();
  void set_translation(const double x, const double y, const double z);
  void set_translation(const double* v);
  void copy(const affine& b);
  void copy_transposed(const affine& b);
  void mult(const affine& b);
  void mult(const affine& b, affine& c);
  void set_rotation(const dMatrix3& rot);
  void set_rotation(const affine& rot);
  void translate(extvec& t);
  void transpose();
  void set_a(const int i, const int j, const double val);
  double get_a(const int i, const int j);
  void mult(extvec& v, extvec& u);
  void subtract(const affine& b);
  double norm();
  void invert_rigidbody();
  void get_translation(extvec& t);
private:
  void print_rows(const int n);
};


class extvec{
  double v[4];
public:
  extvec(){v[3] = 1;}
  extvec(const double x, const double y, const double z);
  inline void set_v(const int i, const double val){v[i] = val;}
  inline double* get_data(){return v;}
  inline double get_v(const int i){return v[i];}
  void set_zeros();
  void print();
  void set(const double x, const double y, const double z);
  void set(const double* a);
  void copy(const extvec& u);
  void copy3(const extvec& u);
  void subtract(const extvec& u);
  double norm();
  void get_components(double& x, double& y, double& z);
  void get_components(double* p);
  void cross(extvec& u, extvec& w);
  void cross(extvec& u);
  void times(const double f);
  void normalize();
  double dot(const extvec& u);
  void to_dvec(dVector3& u);
  void add(const extvec& u);
};

#endif
