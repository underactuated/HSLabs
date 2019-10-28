#include <sstream>
#include <iostream>
#include <math.h>
#include "core.h"

using namespace std;

void str_to_val(const char* str, double* val){
  stringstream ss;
  ss << str;
  while(ss >> *val++){};
}

/*void xmlnode_attr_to_val(xml_node<>* xnode, const char* attr_name, double* val){
  return str_to_val(xnode->first_attribute(attr_name)->value(),val);
  }*/

void xmlnode_attr_to_val(const xml_node<>* xnode, const char* attr_name, double* val){
  xml_attribute<>* attr = xnode->first_attribute(attr_name);
  if(attr){
    str_to_val(xnode->first_attribute(attr_name)->value(),val);
  } else {cout << "WARNING: no attribute " << attr_name << endl;}
}

void xmlnode_attr_to_val(const xml_node<>* xnode, const char* attr_name, string& val){
  xml_attribute<>* attr = xnode->first_attribute(attr_name);
  if(attr){
    val = attr->value();
  } else {cout << "WARNING: no attribute " << attr_name << endl;}
}

double** new_2d_array(int n, int m){
  double** array = new double* [n];
  double** p = array;
  for(int i=0;i<n;i++){
    *p++ = new double [m];
  }
  return array;
}

void delete_2d_array(double** array, int n){
  double** p = array;
  if(p == NULL){return;}
  for(int i=0;i<n;i++){delete [] *p++;}
  delete [] array;
  array = NULL;
}

void save_2d_array(double** array, int n, int m, string fname, bool append_flag){
  ofstream file;
  if (append_flag) {file.open(fname.c_str(),ios_base::app);}
  else {file.open(fname.c_str());}
  double **p = array;
  for(int i=0;i<n;i++){
    double *p1 = *p++;
    for(int j=0;j<m;j++){
      if(j){file << " ";}
      file << *p1++;
    }
    file << endl;
  }
  file.close();
}


double randf(){return double(rand())/RAND_MAX;}


void arrayops::print(double* a){
  for(int i=0;i<n;i++){
    if(i){cout << " ";}
    cout << *a++;
  }
  cout << endl;
}

void arrayops::assign(double* a, const double* a1){
  for(int i=0;i<n;i++){*a++ = *a1++;}
}

double* arrayops::add(double* a, const double* a1){
  double* p = a;
  for(int i=0;i<n;i++){*a++ += *a1++;}
  return p;
}

double* arrayops::subtract(double* a, const double* a1){
  double* p = a;
  for(int i=0;i<n;i++){*a++ -= *a1++;}
  return p;
}

double* arrayops::times(double* a, double b){
  double* p = a;
  for(int i=0;i<n;i++){*a++ *= b;}
  return p;
}

double* arrayops::modulus(double* a, double b){
  double* p = a;
  double bh = b/2;
  for(int i=0;i<n;i++){
    if (*a > bh) {*a -= b;}
    else if (*a < -bh) {*a += b;}
    a++;
  }
  return p;
}

double arrayops::dot(const double* a, const double* a1){
  double s = 0;
  const double*p = a, *p1 = a1;
  for(int i=0;i<n;i++){
    s += (*p++)*(*p1++);
  }
  return s;
}

double arrayops::norm(const double * a){
  return sqrt(dot(a,a));
}

double arrayops::distance(const double* a, const double* a1){
  double* d = new double [n];
  assign(d,a);
  subtract(d,a1);
  double l = norm(d);
  delete [] d;
  return l;
}

void arrayops::assign_scalar(double* a, double b){
  for(int i=0;i<n;i++){*a++ = b;}
}

double arrayops::l1_norm(const double* a){
  double s = 0;
  for(int i=0;i<n;i++){s += fabs(*a++);}
  return s;
}

double** arrayops::new_2d_array(int m){
  return ::new_2d_array(m, n);
}

void arrayops::delete_2d_array(double** array, int m){
  ::delete_2d_array(array, m);
}
