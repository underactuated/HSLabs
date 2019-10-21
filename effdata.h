#include <vector>
#include <map>
#include <list>


struct balltree;

struct efficientdata{
  int data_size;
  vector<double*> points;
  map<double*,int> pis;
  balltree* tree;
  efficientdata(){tree=NULL;}
  void prepare_data(list<vector<double> >& data, int d);
  void prepare_data(double** target_points, int tps_size, int d);
  int get_dim();
  void get_gammat0s(double* a, double b, double c, map<int,double>& gammaa, map<int,double>& t0a);
  void get_tpis(list<int>& tpis, int n, double* btil, double btq0, double btdq0, double sg, double w);
};

typedef pair<double*,double> plane; 

struct balltreetester{
  int dim;
  balltree* tree;
  arrayops ao;
  balltreetester(int dim);
  double* new_rand_point();
  double point_plane_dist(double* pnt, plane pl);
  void test();
  void print_point(double* pnt);
  double norm(double* pnt);
  double dot(double* pnt1, double* pnt2);
};

