#ifndef CORE_H
#define CORE_H

#include "rapidxml-1.13/rapidxml.hpp"
#include "rapidxml-1.13/rapidxml_utils.hpp"

using namespace rapidxml;

void str_to_val(const char* str, double* val);

void xmlnode_attr_to_val(const xml_node<>* xnode, const char* attr_name, double* val);
void xmlnode_attr_to_val(const xml_node<>* xnode, const char* attr_name, std::string& val);

double** new_2d_array(const int n, const int m);
void delete_2d_array(double** array, const int n);
void save_2d_array(double** array, const int n, const int m, const std::string fname, const bool append_flag);

double randf();

// maybe generalize to any type later with templates
class arrayops{
  int n;
public:
  arrayops(const int n_){n = n_;}
  arrayops(){n = 0;}
  void set_n(const int n_){n = n_;}
  void print(double* a);
  void assign(double* a, const double* a1);
  double* add(double* a, const double* a1);
  double* subtract(double* a, const double* a1);
  double* times(double* a, const double b);
  double* modulus(double* a, const double b);
  double dot(const double* a, const double* a1);
  double norm(const double * a);
  double distance(const double* a, const double* a1);
  void assign_scalar(double* a, const double b);
  double l1_norm(const double* a);
  double** new_2d_array(const int m);
  void delete_2d_array(double** array, const int m);
};


#include <iostream>
using namespace std;

template <typename T>
void print_array (const T* a, const int n){
  for(int i=0;i<n;i++){
    if(i){cout << " ";}
    cout << *a++;
  }
  cout << endl;
}

template <typename T>
void print_array (const T* a, const int n, const char* prefix){
  cout << prefix;
  print_array <T> (a,n);
}

#endif
