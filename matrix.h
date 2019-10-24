#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <ode/ode.h>

using namespace std;

class extvec;

class affine{ 
  double a[16]; // matrix is stored column-wise
public:
  void set_zeros();
  void set_unity();
  void print_all();
  void print();
  void set_translation(double x, double y, double z);
  void set_translation(double* v);
  void copy(affine& b);
  void copy_transposed(affine& b);
  void mult(affine& b);
  void mult(affine& b, affine& c);
  inline double* get_data(){return a;}
  void set_rotation(dMatrix3& rot);
  void set_rotation(affine& rot);
  void translate(extvec& t);
  void transpose();
  void set_a(int i, int j, double val);
  double get_a(int i, int j);
  void mult(extvec& v, extvec& u);
  void subtract(affine& b);
  double norm();
  void invert_rigidbody();
  void get_translation(extvec& t);
private:
  void print_rows(int n);
};


class extvec{
  double v[4];
public:
  extvec(){v[3] = 1;}
  extvec(double x, double y, double z);
  void set_zeros();
  void print();
  void set(double x, double y, double z);
  void set(double* a);
  void copy(extvec& u);
  void copy3(extvec& u);
  void subtract(extvec& u);
  double norm();
  inline double* get_data(){return v;}
  void get_components(double& x, double& y, double& z);
  void get_components(double* p);
  void cross(extvec& u, extvec& w);
  void cross(extvec& u);
  void times(double f);
  void normalize();
  double dot(extvec& u);
  void to_dvec(dVector3& u);
  void set_v(int i, double val);
  void add(extvec& u);
  inline double get_v(int i){return v[i];}
};

#endif
