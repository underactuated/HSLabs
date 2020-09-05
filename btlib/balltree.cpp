#include <cstdlib>
#include <cmath>

#include <iostream>

#include "balltree.h"

using namespace std;

void callback::set_dim(int dim_){
  dim=dim_;
  bounds.resize(dim);
}

// computes box enclosing all the points in psv
void callback::compute_bounds(vector<double*>& psv){
  //if(psv.size()==0){cout<<"error, psv size 0"<<endl;exit(1);}
  for(int i=0;i<dim;i++){
    vector<double*>::iterator it = psv.begin();
    double lb=1e10, ub=-1e10;
    for(;it!=psv.end();it++){
      double x = (*it)[i];
      if(x<lb){lb=x;}
      if(x>ub){ub=x;}
    }
    bounds[i]=make_pair(lb,ub);
  }
}

// sets r as half diagonal of the enclosing box
void callback::set_xcr(double*& xc, double& r){
  double s=0;
  for(int i=0;i<dim;i++){
    double lb = bounds[i].first, ub = bounds[i].second;
    xc[i]=(lb+ub)/2.;
    double delb = ub-lb;
    s+=delb*delb;
  }
  r=sqrt(s)/2.;
}

// index of the largest enclosing box's side
void callback::idim_maxspread(){
  double delbmax=0;
  for(int i=0;i<dim;i++){
    double lb = bounds[i].first, ub = bounds[i].second;
    double delb = ub-lb;
    if(delb>delbmax){delbmax=delb; idim=i;}
  }
}

void callback::print(){
  cout << "dim = " << dim << endl;
  cout << "bounds:";
  for(int i=0;i<dim;i++){
    cout << " (" << bounds[i].first << "," << bounds[i].second << ")";
  }
  cout << endl;
  cout << "idim = " << idim << endl;
}

/*// SEEMS INCORRECT
// should be lb = (dotprod(a,n->xc,dim)+b)/norma-n->r, etc
pair<float,float> callback::planedist_lubounds(bnode* n){
  float lb = dotprod(a,n->xc,dim)+b-norma*n->r;
  float ub = dotprod(a,n->xc,dim)+b+norma*n->r;
  if (lb<0&&ub>0) {
    if (lb+ub>0) {lb=0;} else {ub=-lb; lb=0;}
  } else if (ub<0) {float tmp=ub; ub=-lb; lb=-tmp;}
  return make_pair(lb,ub);
  }*/

// CORRECTED VERSION 
// min and max ball-plane distance
pair<float,float> callback::planedist_lubounds(bnode* n){
  float lb = (dotprod(a,n->xc,dim)+b)/norma-n->r;
  float ub = (dotprod(a,n->xc,dim)+b)/norma+n->r;
  if (lb<0&&ub>0) {
    if (lb+ub>0) {lb=0;} else {ub=-lb; lb=0;}
  } else if (ub<0) {float tmp=ub; ub=-lb; lb=-tmp;}
  return make_pair(lb,ub);
}

// sets plane params
void callback::set_ab(double* a_, double b_){
  a=a_;
  norma = norm(a,dim);
  b=b_;
}

//int jj=0;
void callback::knnplane(bnode* n){
  //  jj++;
  if(n->n){
    pair<float,float> lubs = planedist_lubounds(n);
    if(lubs.first<=klps->bound||klps->if_mlessk()){
      klps->add_node(n,lubs.first,lubs.second);
    }
    return;
  }
  
  bnode* n0 = n->left_child;
  bnode* n1 = n->right_child;
  pair<float,float> lubs = planedist_lubounds(n0);
  float lb0 = lubs.first;
  lubs = planedist_lubounds(n1);
  float lb1 = lubs.first;
  if(lb1<lb0){n1=n0; n0=n->right_child; lb1=lb0;}
  knnplane(n0);
  //cout<<klps->m<<" "<<klps->bound<<" "<<lb1<<endl;
  if(klps->if_mlessk()){knnplane(n1);}
  else if(lb1<=klps->bound){knnplane(n1);}
}

void callback::klps_nns(map<double,double*>& nns){
  list<bnode*> nodes;
  klps->get_nodes(nodes);
  map<double,double*> nns_pos, nns_neg;
  list<bnode*>::iterator it = nodes.begin();
  for(;it!=nodes.end();it++){
    double** ps = (*it)->ps;
    int n = (*it)->n;
    for(int i=0;i<n;i++){
      double* p = *ps++;
      double d = dotprod(a,p,dim)+b;
      if(d>0){nns_pos[d]=p;}
      else{nns_neg[-d]=p;}
    }
  }
  int k = klps->k;
  int posneg_size = nns_pos.size()+nns_neg.size();
  if(posneg_size<k){k=posneg_size;}
  nns_pos[1e10]=NULL;
  nns_neg[1e10]=NULL;
  map<double,double*>::iterator it1 = nns_pos.begin(), it2 = nns_neg.begin();
  for(int i=0;i<k;i++){
    if((*it1).first<(*it2).first){
      nns[(*it1).first]=(*it1).second;
      it1++;
    } else {
      nns[-(*it2).first]=(*it2).second;
      it2++;
    }
  }
  //cout<<nns.size()<<endl;
  //if(nns.size()<k){cout<<"size < k"<<endl;exit(1);}
}


bnode::bnode(callback& cb){
  xc = new double [cb.dim];
  n=0;
}

bnode::~bnode(){
  delete [] xc;
  if(n){delete [] ps;}
  else {
    delete left_child;
    delete right_child;
  }
}

void bnode::add_points(vector<double*>& psv, callback cb){
  cb.compute_bounds(psv);
  cb.set_xcr(xc,r);
  if (psv.size() > MAX_BLEAF_SIZE){
    make_children(psv,cb);
  } else {
    n=psv.size();
    ps = new double* [n];
    for(int i=0;i<n;i++){ps[i]=psv[i];}
  }
}

void bnode::make_children(vector<double*>& psv, callback cb){
  cb.idim_maxspread();
  int imax = cb.idim;
  vector<double*> left_psv, right_psv;
  vector<double*>::iterator it = psv.begin();
  for(;it!=psv.end();it++){
    if((*it)[imax]<xc[imax]){left_psv.push_back(*it);}
    else{right_psv.push_back(*it);}
  }
  if(left_psv.size()==0){
    left_psv.push_back(right_psv.back());
    right_psv.pop_back();
  }
  left_child = new bnode (cb);
  right_child = new bnode (cb);
  left_child->add_points(left_psv,cb);
  right_child->add_points(right_psv,cb);
  //if(n){cout<<"error"<<endl;exit(1);} // temp
}

void bnode::print(callback& cb){
  cout << "xc: ";
  print_array(xc,cb.dim);
  cout << "r = " << r << endl;
  cout << "n = " << n << endl;
}

pair<int,double**> bnode::get_points(){
  return make_pair<int,double**> (n,ps);
}

bnode* bnode::get_child(int i){
  if(i){return right_child;}
  else{return left_child;}
}


balltree::balltree(int dim){
  cb.set_dim(dim);
  root = new bnode (cb);
}

balltree::~balltree(){delete root;}

void balltree::add_points(vector<double*>& psv){
  cb.compute_bounds(psv);
  cb.idim_maxspread();
  root->add_points(psv,cb);
}

void balltree::knnplane(double* a, double b, int k, map<double,double*>& nns){
  cb.set_ab(a,b);
  klowestpoints klps(k);
  cb.klps=&klps;
  cb.knnplane(root);
  //klps.print();cout<<"jj="<<jj<<endl;
  cb.klps_nns(nns);
}

void balltree::knnplane1(double* a, double b, int k, int m, map<double,double*>& nns){
  goaldistance dist (get_dim());
  dist.set(a,b);
  approxsearch as;
  as.knn(root,k,&dist,k+m,nns);
}

void balltree::print(){
  cout << "balltree:" << endl;
  cb.print();
}


// notice: for m<=k, m1==m; for m>k, k<=m1<=m.
void klowestpoints::add_node(bnode* n, float lbound, float ubound){
  if(m<k){
    lbounds.push(lbound,n);
    m+=n->n;
    ubounds.push(ubound,n);
    m1+=n->n;
  } else if (lbound<=bound){
    lbounds.push(lbound,n);
    m+=n->n;
    if(ubound<=bound){
      ubounds.push(ubound,n);
      m1+=n->n;
    }
  }
  adjust_bounds();
}

void klowestpoints::adjust_bounds(){
  if(m1>k){
    int n1 = ubounds.val()->n;
    if(m1-n1 >= k){
      m1-=n1;
      ubounds.pop();
    }
  }
  bound = ubounds.key();
  while(lbounds.key()>bound){
    m-=lbounds.val()->n;
    lbounds.pop();
  }
}

void klowestpoints::print(){
  cout << "k = " << k <<" m = " << m << " m1 = " << m1 << endl;
  cout << "bound = " << bound << endl;
  cout << "lbounds: " << lbounds.size() << endl;
  cout << "ubounds: " << ubounds.size() << endl;
}

void klowestpoints::get_nodes(list<bnode*>& nodes){
  while(lbounds.size()){
    nodes.push_back(lbounds.val());
    lbounds.pop();
  }
} 


double frand(){return double(rand()%RAND_MAX)/RAND_MAX;}

double* random_array(int dim){
  double* p = new double [dim];
  for(int i=0;i<dim;i++){p[i]=frand()*2-1;}
  return p;
}

void tree_stat(bnode* n){
  int stat[]={0,0,0};
  string fields[]={"points   ","int nodes","leaves   "};
  collect_stat(n,stat);
  cout << "---------------" << endl;
  for(int i=0;i<3;i++){
    cout << fields[i] << " = " << stat[i] << endl;
  }
}

void collect_stat(bnode* n, int* stat){
  if(n->n>0){
    stat[0]+=n->n;
    stat[2]++;
  } else {
    stat[1]++;
    collect_stat(n->left_child,stat);
    collect_stat(n->right_child,stat);
  }
}

void all_leaves(bnode* n, list<bnode*>& nodes){
  if(n->n>0){nodes.push_back(n);}
  else {
    all_leaves(n->left_child,nodes);
    all_leaves(n->right_child,nodes);
  }
}

void print_array(double* a, int dim){
  cout << *a++;
  for(int i=1;i<dim;i++){cout << " " << *a++;}
  cout << endl;
}

/*double dotprod(double* a, double* b, int dim){
  double s=0;
  for(int i=0;i<dim;i++){s+=(*a++)*(*b++);}
  return s;
  }*/

double dotprod(double* a, double* b, int dim){
  double s=0;
  double *pa = a, *pb = b;
  for(int i=0;i<dim;i++){s+=(*pa++)*(*pb++);}
  return s;
}

double norm(double* a, int dim){
  return sqrt(dotprod(a,a,dim));
}

void normalize(double* a, int dim){
  double z = norm(a,dim);
  for(int i=0;i<dim;i++){*a++/=z;}
}


void approxsearch::knn (bnode* root, int k_, goaldistance* dist_, int m_, map<double,double*>& nns){
  k = k_;
  m = m_;
  if (m < k) {cout<<"ERROR: m < k"<<endl;exit(1);}
  dist = dist_;
  m1 = 0;
  max_dist = 1e10;
  score_node_heap.clear();
  dist_point_heap.clear();
  process_bnode(root);
  //cout<<"heap size="<<dist_point_heap.size()<<endl;
  /*for(int i=0;i<k;i++){
    double d = dist_point_heap.key();
    nns[d] = dist_point_heap.pull_val();
    }*/
  while(!dist_point_heap.is_empty()){
    double d = dist_point_heap.key();
    nns[d] = dist_point_heap.pull_val();
  }
  //cout<<"k="<<k<<" m="<<m<<" m1="<<m1<<endl;
}

void approxsearch::process_bnode(bnode* node){
  //if children, score them, pick best from score_node_heap, process it
  // else check points
  if(node->is_leaf()){
    check_points(node->get_points());
  } else {
    for(int i=0;i<2;i++){
      bnode* child = node->get_child(i);
      to_score_node_heap(child);
    }
  }
  while((m1 < m) && !score_node_heap.is_empty()){
    bnlb elem = score_node_heap.pull_val();
    if(elem.second < max_dist){
      process_bnode(elem.first);
    }
  }
}

void approxsearch::to_score_node_heap(bnode* node){
  double* xc = node->get_xc();
  double r = node->get_r();
  double d = dist->dist(xc);
  double ub = d+r;
  double lb = (r < d)? d-r : 0;
  double score = (ub+lb)/2;
  bnlb elem (node,lb);
  score_node_heap.insert(-score,elem);
}

void approxsearch::check_points(pair<int,double**> points){
  int n = points.first;
  double** ps = points.second;
  for(int i=0;i<n;i++){
    double* p = *ps++;
    double d = dist->dist(p);
    if(d < max_dist){
      dist_point_heap.insert(d,p);
      if(dist_point_heap.size() > k){
	dist_point_heap.pull_val();
	max_dist = dist_point_heap.key();
      }
      //cout<<"max_dist="<<max_dist<<endl;
    }
    m1++;
  }
  //cout<<"m1="<<m1<<endl;
}

/*
double approxsearch::score_bnode(bnode* node){
  double* xc = node->get_xc();
  double r = node->get_r();
  double d = dist->dist(xc);
  if(r < d){return d;}
  else {return (d+r)/2;}
}
*/

double goaldistance::dist(double* x){
  return fabs(dotprod(x,a,dim)+b)/norm(a,dim);
}

void goaldistance::set_random(){
  a = random_array(dim);
  b = frand()*2-1;
}


void random_btree(int n, balltree& tree){
  vector<double*> psv;
  int dim = tree.get_dim();
  for(int i=0;i<n;i++){
    double* p = random_array(dim);
    psv.push_back(p);
  }
  tree.add_points(psv);
}
