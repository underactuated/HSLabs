/////////////////////////////////////////////////
// effdata.h: declaration of classes implementing
// approximate search for efficient search of
// target points in by CPC controller.
/////////////////////////////////////////////////
#ifndef EFFDATA_H
#define EFFDATA_H

#include <map>
#include <list>

#include "core.h"

class balltree;

class efficientdata{
  int data_size;
  vector<double*> points;
  map<double*,int> pis;
  balltree* tree;
public:
  efficientdata(){tree=NULL;}
  void prepare_data(const list<vector<double> >& data, int d);
  void prepare_data(double** target_points, int tps_size, int d);
  int get_dim();
  void get_gammat0s(double* a, double b, double c, map<int,double>& gammaa, map<int,double>& t0a);
  void get_tpis(list<int>& tpis, int n, double* btil, double btq0, double btdq0, double sg, double w);
};

typedef pair<double*,double> plane; 

class balltreetester{
  int dim;
  balltree* tree;
  arrayops ao;
public:
  balltreetester(int dim);
  double* new_rand_point();
  double point_plane_dist(const double* pnt, plane pl);
  void test();
  void print_point(double* pnt);
  double norm(const double* pnt);
  double dot(const double* pnt1, const double* pnt2);
};

#endif
