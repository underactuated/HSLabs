#include <math.h>
#include "core.h"
#include "effdata.h"
#include "balltreelib1/balltree.h"


void efficientdata::prepare_data(list<vector<double> >& data, int d){
  data_size=0;
  list<vector<double> >::iterator it = data.begin();
  for(;it!=data.end();it++){
    if((int)points.size()<=data_size){
      points.push_back(new double [2*d]);
      pis[points.back()]=data_size;
    }
    double* p = points[data_size];
    vector<double>::iterator it1 = (*it).begin();
    for(int i=0;i<2*d;i++){*p++ = (*it1++);}
    data_size++;
  }
  if((int)points.size()>data_size){cout<<"WARNING: points size > data size: "<<points.size()<<" "<<data_size<<endl;}
  if(tree){delete tree;}
  tree = new balltree (d);
  tree->add_points(points);
  //tree_stat(tree->root);exit(1);
  //tree->print();exit(1);
}

void efficientdata::prepare_data(double** target_points, int tps_size, int d){
  data_size=0;
  points.clear();
  pis.clear();
  double **p = target_points;
  for(int i=0;i<tps_size;i++){
    for(int j=0;j<2;j++){(*p)[j] = 0;} // experim
    points.push_back(*p++);
    pis[points.back()] = data_size++;
  }
  if(tree){delete tree;}
  tree = new balltree (d);
  tree->add_points(points);
  //tree_stat(tree->root);exit(1);
  //tree->print();exit(1);
}

int efficientdata::get_dim(){return tree->cb.dim;}

void efficientdata::get_gammat0s(double* a, double b, double c, map<int,double>& gammaa, map<int,double>& t0a){
  int k = 1000/2;
  int d = get_dim();
  map<double,double*> nns;
  //tree->knnplane(a,b,k,nns);
  tree->knnplane1(a,b,k,2*k+100,nns);
  map<double,double*>::iterator it = nns.begin();
  for(;it!=nns.end();it++){
    double* p = (*it).second;
    int i = pis[p];
    gammaa[i] = c/dotprod(a,p+d,d);
    t0a[i] = (*it).first/c;
  }
}

void efficientdata::get_tpis(list<int>& tpis, int n, double* btil, double btq0, double btdq0, double sg, double w){
  map<int,double> gammaa, t0a;
  get_gammat0s(btil,-btq0,btdq0,gammaa,t0a);
  
  map<double,int> loss_tpi;
  map<int,double>::iterator it = t0a.begin();
  for(;it!=t0a.end();it++){
    int tpi = (*it).first;
    double t0 = (*it).second;
    double s = 1./gammaa[tpi];
    double a = w*t0, a1 = s-sg; 
    double loss = a*a + a1*a1;
    loss_tpi[loss] = tpi;
  }

  map<double,int>::iterator it1 = loss_tpi.begin();
  int size = loss_tpi.size();
  if(n > size){n = size;}
  for(int i=0;i<n;i++){
    tpis.push_back((*it1).second);
    //cout << tpis.back() << endl;
    it1++;
  }
}


balltreetester::balltreetester(int dim_){
  dim = dim_;
  tree = new balltree (dim);
  ao.set_n(dim);
}

void balltreetester::test(){
  double* pnt = new_rand_point();
  print_point(pnt);
  plane pl (pnt,randf());
  double* pnt1 = new_rand_point();
  cout << point_plane_dist(pnt1,pl) << endl;
  cout << "b=" << pl.second << endl;

  vector<double*> points;
  for(int i=0;i<1000;i++){
    double* pnt = new_rand_point();
    points.push_back(pnt);
  }

  tree->add_points(points);

  int k = 10;
  map<double,double*> nns;
  //tree->knnplane(pl.first,pl.second,k,nns);
  tree->knnplane1(pl.first,pl.second,k,200,nns);
  map<double,double*>::iterator it = nns.begin();
  for(;it!=nns.end();it++){
    double* p = (*it).second;
    cout << (*it).first << " " << point_plane_dist(p,pl) << endl;
  }
}

double* balltreetester::new_rand_point(){
  double* pnt = new double [dim];
  double* p = pnt;
  for(int i=0;i<dim;i++){*p++ = randf()*2.-1;}
  return pnt;
}

void balltreetester::print_point(double* pnt){
  print_array<double>(pnt,dim);
}

double balltreetester::point_plane_dist(double* pnt, plane pl){
  double* a = pl.first;
  double b = pl.second;
  return fabs(dot(pnt,a)+b)/norm(a);
}

double balltreetester::norm(double* pnt){
  return ao.norm(pnt);
}

double balltreetester::dot(double* pnt1, double* pnt2){
  return ao.dot(pnt1,pnt2);
}
