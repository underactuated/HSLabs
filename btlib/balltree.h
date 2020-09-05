#include <list>
#include <vector>
#include <map>
#include <algorithm>

#include "heap.h"

namespace myheap{
#include "myheap.h"
}

const unsigned int MAX_BLEAF_SIZE=2;//5;

using namespace std;

struct bnode;
struct klowestpoints;

struct callback{
  int dim;
  vector<pair<double,double> > bounds;
  int idim;
  double *a, norma, b;
  klowestpoints *klps;
  void set_dim(int dim_);
  void compute_bounds(vector<double*>& ps);
  void set_xcr(double*& xc, double& r);
  void idim_maxspread();
  void print();
  void set_ab(double* a_, double b_);
  pair<float,float> planedist_lubounds(bnode* n);
  void knnplane(bnode* n);
  //bool if_leftcloser(bnode* n);
  void klps_nns(map<double,double*>& nns);
  int get_dim(){return dim;}
};

struct bnode{
  double *xc, r;
  bnode *left_child, *right_child;
  int n;
  double** ps;
  bnode(callback& cb);
  ~bnode();
  void add_points(vector<double*>& psv, callback cb);
  void make_children(vector<double*>& psv, callback cb);
  void print(callback& cb);
  bool is_leaf(){return (n>0);}
  pair<int,double**> get_points();
  bnode* get_child(int i);
  double* get_xc(){return xc;}
  double get_r(){return r;}
};

struct balltree{
  callback cb;
  bnode* root;
  balltree(int dim);
  ~balltree();
  void add_points(vector<double*>& psv);
  void print_node(bnode* n){n->print(cb);}
  void knnplane(double* a, double b, int k, map<double,double*>& nns);
  void knnplane1(double* a, double b, int k, int m, map<double,double*>& nns);
  void print();
  bnode* get_root(){return root;}
  int get_dim(){return cb.get_dim();}
};

typedef myheap::heap<bnode> heapbn_t;

struct klowestpoints{
  int k, m, m1; // m all points, at least m1 points below bound
  float bound;
  heapbn_t lbounds, ubounds;
  klowestpoints(int k_){k=k_; m=0; m1=0; bound=0;}
  void add_node(bnode* n, float lbound, float ubound);
  void adjust_bounds();
  void print();
  bool if_mlessk(){return (m<k);}
  void get_nodes(list<bnode*>& nodes);
};

double* random_array(int dim);
void tree_stat(bnode* n);
void collect_stat(bnode* n, int* stat);
void all_leaves(bnode* n, list<bnode*>& nodes);
void print_array(double* a, int dim);
double dotprod(double* a, double* b, int dim);
double norm(double* a, int dim);
void normalize(double* a, int dim);
double frand();

struct goaldistance;

typedef pair<bnode*,double> bnlb; // bnlb = bnode lower bound

struct approxsearch{
  //heap<double,bnode*> score_node_heap;
  heap<double,bnlb> score_node_heap;
  heap<double,double*> dist_point_heap;
  int k, m, m1; // m points to check, m1 points checked
  goaldistance* dist;
  double max_dist;
  void knn (bnode* root, int k, goaldistance* dist, int m, map<double,double*>& nns);
  void process_bnode(bnode* node);
  void to_score_node_heap(bnode* node);
  void check_points(pair<int,double**> points);
  //double score_bnode(bnode* node);
};

struct goaldistance{
  int dim;
  double *a, b;
  goaldistance(int dim_){dim = dim_;}
  double dist(double* x);
  void set_random();
  void set(double* a_, double b_){a = a_; b = b_;}
};

void random_btree(int n, balltree& tree);



