/////////////////////////////////////////////////
// core.h: declaration of some useful functions
// and a class for array manipulation
/////////////////////////////////////////////////
#ifndef CORE_H
#define CORE_H

#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_utils.hpp"

using namespace rapidxml;

void str_to_val(const char* str, double* val);

void xmlnode_attr_to_val(const xml_node<>* xnode, const char* attr_name, double* val);
void xmlnode_attr_to_val(const xml_node<>* xnode, const char* attr_name, std::string& val);

double** new_2d_array(int n, int m);
void delete_2d_array(double** array, int n);
void save_2d_array(double** array, int n, int m, std::string fname, bool append_flag);

double randf();

void save_vector_of_records(std::vector<double*>& recs, int rec_len, int offset, std::string fname);
void load_vector_of_records(std::vector<double*>& recs, int rec_len, std::string fname);

// Class arrayops ( = array operations) object simplifis certain array 
// operations on arrays of given length n, such as vector arifmetic.
// TODO: maybe generalize to any type with templates?
class arrayops{
  int n;  // array size
  double* tmp; // extra array
public:
  arrayops(int n);
  arrayops();
  ~arrayops();
  void set_n(int n_);
  void print(double* a);
  void assign(double* a, const double* a1) const;
  double* add(double* a, const double* a1) const;
  double* subtract(double* a, const double* a1) const;
  double* times(double* a, double b) const;
  double* modulus(double* a, double b) const;
  double dot(const double* a, const double* a1) const;
  double norm(const double * a) const;
  double distance(const double* a, const double* a1);
  void assign_scalar(double* a, double b) const;
  double l1_norm(const double* a) const;
  double** new_2d_array(int m) const;
  void delete_2d_array(double** array, int m) const;
};


#include <iostream>
using namespace std;

template <typename T>
void print_array (const T* a, int n){
  for(int i=0;i<n;i++){
    if(i){cout << " ";}
    cout << *a++;
  }
  cout << endl;
}

template <typename T>
void print_array (const T* a, int n, const char* prefix){
  cout << prefix;
  print_array <T> (a,n);
}

template <typename T>
void print_vector (const vector<T>& v) {
  typename vector<T>::const_iterator it = v.begin();
  for(;it!=v.end();it++){
    if(it != v.begin()){cout << " ";}
    cout << *it;
  }
  cout << endl;
}

template <typename T>
void print_vector (const vector<T>& v, const char* prefix) {
  cout << prefix;
  print_vector <T> (v);
}

template <typename T>
void load_array (T* a, int n, string fname){
  ifstream file;
  file.open(fname.c_str());
  for(int i=0;i<n;i++){file >> *a++;}
  file.close();
}


#endif
