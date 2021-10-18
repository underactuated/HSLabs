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

// Allocates  2d (two-dimensional, n x m) array.
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

// Saves 2d array to file fname, appending if append_flag.
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

// uniform sampling from [0,1]
double randf(){return double(rand())/RAND_MAX;}

// record is an array of double
void save_vector_of_records(vector<double*>& recs, int rec_len, int offset, string fname){
  int size = recs.size();
  double** recs_array = new double* [size];
  std::copy(recs.begin(),recs.end(),recs_array);
  save_2d_array(recs_array+offset,size-offset,rec_len,fname,false);
  delete [] recs_array;
}

void load_vector_of_records(vector<double*>& recs, int rec_len, string fname){
  ifstream file;
  file.open(fname.c_str());
  string str;
  while(getline(file,str)){
    double* rec = new double [rec_len];
    stringstream ss; ss << str;
    for(int i=0;i<rec_len;i++){
      if(!(ss >> rec[i])){
	cout<<"ERROR: str too short in load_vector_of_records"<<endl;exit(1);
      }
    }
    recs.push_back(rec);
  }
}


arrayops::arrayops(int n_){
  tmp = NULL;
  set_n(n_);
}

arrayops::arrayops(){
  tmp = NULL;
  n = 0;
}

arrayops::~arrayops(){
  if(tmp){delete [] tmp;}
}

void arrayops::set_n(int n_){
  n = n_;
  if(tmp){delete [] tmp; tmp = NULL;}
}

void arrayops::print(double* a){
  for(int i=0;i<n;i++){
    if(i){cout << " ";}
    cout << *a++;
  }
  cout << endl;
}

// a = a1
void arrayops::assign(double* a, const double* a1) const {
  for(int i=0;i<n;i++){*a++ = *a1++;}
}

// a += a1
double* arrayops::add(double* a, const double* a1) const {
  double* p = a;
  for(int i=0;i<n;i++){*a++ += *a1++;}
  return p;
}

// a -= a1
double* arrayops::subtract(double* a, const double* a1) const {
  double* p = a;
  for(int i=0;i<n;i++){*a++ -= *a1++;}
  return p;
}

// a *= b, where b is scalar
double* arrayops::times(double* a, double b) const {
  double* p = a;
  for(int i=0;i<n;i++){*a++ *= b;}
  return p;
}

// modulus is used for arguments of periodic functions (with period b)
// to bring a to (-b/2,b/2]. It is assumed that a is in (-3*b/2,3*b/2].
double* arrayops::modulus(double* a, double b) const {
  double* p = a;
  double bh = b/2;
  for(int i=0;i<n;i++){
    if (*a > bh) {*a -= b;}
    else if (*a <= -bh) {*a += b;}
    a++;
  }
  return p;
}

// dot-product of a and a1
double arrayops::dot(const double* a, const double* a1) const {
  double s = 0;
  const double*p = a, *p1 = a1;
  for(int i=0;i<n;i++){
    s += (*p++)*(*p1++);
  }
  return s;
}

// l2 norm
double arrayops::norm(const double * a) const {
  return sqrt(dot(a,a));
}

// norm(a-a1)
double arrayops::distance(const double* a, const double* a1){
  if(tmp == NULL){tmp = new double [n];}
  assign(tmp,a);
  subtract(tmp,a1);
  double l = norm(tmp);
  return l;
}

// a = b, where b is scalar
void arrayops::assign_scalar(double* a, double b) const {
  for(int i=0;i<n;i++){*a++ = b;}
}

double arrayops::l1_norm(const double* a) const {
  double s = 0;
  for(int i=0;i<n;i++){s += fabs(*a++);}
  return s;
}

double** arrayops::new_2d_array(int m) const {
  return ::new_2d_array(m, n);
}

void arrayops::delete_2d_array(double** array, int m) const {
  ::delete_2d_array(array, m);
}
