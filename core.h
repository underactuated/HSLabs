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
  void assign(double* a, const double* a1);
  double* add(double* a, const double* a1);
  double* subtract(double* a, const double* a1);
  double* times(double* a, double b);
  double* modulus(double* a, double b);
  double dot(const double* a, const double* a1);
  double norm(const double * a);
  double distance(const double* a, const double* a1);
  void assign_scalar(double* a, double b);
  double l1_norm(const double* a);
  double** new_2d_array(int m);
  void delete_2d_array(double** array, int m);
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

#endif
