#include <math.h>
#include "matrix.h"


void affine::set_zeros(){
  double *p = a;
  for(int i=0;i<16;i++){*p++ = 0;}
}

void affine::set_unity(){
  set_zeros();
  double *p = a;
  for(int i=0;i<4;i++){
    *p = 1;
    p += 5;
  }
}

void affine::set_translation(double x, double y, double z){
  set_unity();
  double *p = a + 12;
  *p++ = x;
  *p++ = y;
  *p++ = z;  
}

void affine::set_translation(const double* v){
  return set_translation(*v,*(v+1),*(v+2));
}

void affine::print_rows(int n) const {
  affine b;
  b.copy_transposed(*this);
  double *p = b.a;
  for(int i=0;i<n;i++){
    i==0 ? cout << "[" : cout << " ";
    cout << *p++;
    for(int j=0;j<3;j++){cout << " " << *p++;}
    if (i==n-1) {cout << "]";}
    cout << endl;
  }
}

// print all rows for debugging
void affine::print_all() const {
  print_rows(4);
}

// prints rotation and translation (first 3 rows)
void affine::print() const {
  print_rows(3);
}

// copies from b
void affine::copy(const affine& b){
  double *p = a;
  const double *p1 = b.a;
  for(int i=0;i<16;i++){*p++ = *p1++;}
}

void affine::copy_transposed(const affine& b){
  double *p = a;
  const double *p1 = b.a;
  for(int i=0;i<4;i++){
    p1 = b.a + i;
    for(int j=0;j<4;j++){
      *p++ = *p1;
      p1+=4;
    }
  }
}

void affine::mult(const affine& b){
  affine c;
  c.copy_transposed(*this);
  double *p = a, *pc;
  const double *pb = b.a;
  for(int i=0;i<4;i++){
    pc = c.a;
    for(int j=0;j<4;j++){
      double s = 0;
      double *p1 = pc;
      const double *p2 = pb;
      for(int k=0;k<4;k++){
	s += (*p1++)*(*p2++);
      }
      *p++ = s;
      pc += 4;
    }
    pb +=4;
  }
}

void affine::mult(const affine& b, affine& c) const {
  c.copy(*this);
  c.mult(b);
}

void affine::set_rotation(const dMatrix3& rot){
  double* p = a;
  const dReal* p1 = rot;
  for(int i=0;i<12;i++){*p++ = *p1++;}
  for(int i=0;i<3;i++){*p++ = 0;}
  *p = 1;
}

void affine::set_rotation(const affine& rot){
  double *p = a;
  const double *p1 = rot.a;
  for(int i=0;i<12;i++){*p++ = *p1++;}
  for(int i=0;i<3;i++){*p++ = 0;}
  *p = 1;
}

void affine::translate(const extvec& t){
  double *p = a+12;
  const double *p1 = t.get_data();
  for(int i=0;i<3;i++){*p++ += *p1++;}
}

void affine::transpose(){
  affine A;
  A.copy_transposed(*this);
  copy(A);
}

// i and j are counted from 0
void affine::set_a(int i, int j, double val){
  double* p = a + (j*4+i);
  *p = val;
}
 
double affine::get_a(int i, int j){
  return *(a + (j*4+i));
}

void affine::mult(const extvec& v, extvec& u) const {
  //double *p, *p2 = u.get_data();
  double *p2 = u.get_data();
  for(int i=0;i<4;i++){
    //p = a+i;
    const double *p = a+i;
    const double *p1 = v.get_data();
    double s = 0;
    for(int j=0;j<4;j++){
      s += (*p)*(*p1++);
      p += 4;
    }
    *p2++ = s;
  }
}

void affine::subtract(const affine& b){
  double *p = a;
  const double *p1 = b.a;
  for(int i=0;i<15;i++){*p++ -= *p1++;}
}

// Frobenius norm
double affine::norm() const {
  const double *p = a;
  double s = 0;
  for(int i=0;i<15;i++){s += (*p)*(*p); p++;}
  return sqrt(s);
}

void affine::invert_rigidbody(){
  extvec t, t1;
  double *p = a + 12, *p1 = t.get_data();
  for(int i=0;i<3;i++){
    *p1++ = -(*p);
    *p++ = 0;
  }
  transpose();
  mult(t,t1);
  translate(t1);
}
 
void affine::get_translation(extvec& t) const {
  double *p = t.get_data();
  const double *p1 = a + 12;
  for(int i=0;i<3;i++){*p++ = *p1++;}
}


extvec::extvec(double x, double y, double z){
  v[3] = 1;
  this->set(x,y,z);
}
 
void extvec::set_zeros(){
  double* p = v;
  for(int i=0;i<3;i++){*p++ = 0;}
}

void extvec::print() const {
  const double* p = v;
  for(int i=0;i<3;i++){
    (i == 0)? cout << "[" : cout << " ";
    cout << *p++;
  }
  cout << "]" << endl;
}
 
void extvec::set(double x, double y, double z){
  double* p = v;
  *p++ = x;
  *p++ = y;
  *p = z;
}
 
void extvec::set(const double* a){
  double* p = v;
  for(int i=0;i<3;i++){*p++ = *a++;}
}
 
void extvec::copy(const extvec& u){
  double *p = v;
  const double *p1 = u.v;
  for(int i=0;i<4;i++){*p++ = *p1++;}
}
 
void extvec::copy3(const extvec& u){
  set(u.v);
}
 
void extvec::subtract(const extvec& u){
  double *p = v;
  const double *p1 = u.v;
  for(int i=0;i<4;i++){*p++ -= *p1++;}
}

double extvec::norm() const {
  const double *p = v;
  double s = 0;
  for(int i=0;i<3;i++){s += (*p)*(*p); p++;}
  return sqrt(s);
}

void extvec::get_components(double& x, double& y, double& z) const {
  const double *p = v;
  x = *p++;
  y = *p++;
  z = *p;
}

void extvec::get_components(double* p) const {
  const double *p0 = v;
  for(int i=0;i<3;i++){*p++ = *p0++;}
}

void extvec::cross(const extvec& u, extvec& w) const {
  double x1, y1, z1, x2, y2, z2;
  get_components(x1,y1,z1);
  u.get_components(x2,y2,z2);
  w.set(y1*z2-z1*y2,z1*x2-x1*z2,x1*y2-y1*x2);
}

void extvec::cross(const extvec& u){
  extvec w;
  cross(u,w);
  copy3(w);
}

void extvec::times(double f){
  double *p = v;
  for(int i=0;i<3;i++){*p++ *= f;}
}

void extvec::normalize(){
  double* p = v;
  double len = norm();
  for(int i=0;i<3;i++){*p++ /= len;}
}

double extvec::dot(const extvec& u) const {
  const double *p = v, *p1 = u.v;
  double s = 0;
  for(int i=0;i<3;i++){s += (*p++)*(*p1++);}
  return s;
}

void extvec::to_dvec(dVector3& u) const {
  dReal *p = u; 
  const double *p1 = v;
  for(int i=0;i<3;i++){*p++ = *p1++;}
}

void extvec::add(const extvec& u){
  double *p = v;
  const double *p1 = u.v;
  for(int i=0;i<3;i++){*p++ += *p1++;}
}
