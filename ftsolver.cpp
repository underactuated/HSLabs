//#include <map>
//#include "core.h"
//#include "matrix.h"
//#include "visualization.h"
//#include "model.h"
//#include "lik.h"
//#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include "pergen.h"
//#include "periodic.h"
//#include "dynrec.h"
#include "ftsolver.h"

using namespace Eigen;

forcetorquesolver::forcetorquesolver(const periodic* per){
  n = per->get_number_of_dynparts();
  nf = per->get_nfeet();
  parentis = per->get_parentis();
  footis = per->get_footis();
  masses = per->get_masses();
  mask_l0 = 0;
  jzaxis = NULL;
  per->hinge_joint_part_ids(jpart_ids);
}

forcetorquesolver::~forcetorquesolver(){
  if(!jzaxis){delete [] jzaxis;}
}

void forcetorquesolver::solve_forcetorques(dynrecord* dynrec){
  VectorXd x, y;
  solve_forcetorques(dynrec,x,y);

  cout << "x:" <<endl;
  cout << x.head(15) << endl;
}

void  check_contact_forces(const VectorXd& x0, const MatrixXd& N, const VectorXd& c, double f, const VectorXd& y1, list<pair<int,int> >& mask){

  VectorXd d = c;
  d = d.array().square();

  list<pair<int,int> >::iterator it = mask.begin();
  for(;it!=mask.end();it++){
    int i0 = (*it).first, imax = (*it).second;
    for(int i=i0;i<imax;i++){d(i) *= f;}
  }

  //cout<<d.transpose()<<endl;exit(1);
  DiagonalMatrix<double, Dynamic> D (d);
  MatrixXd m = N.transpose()*D*N;
  //cout<<"det = "<<m.determinant()<<endl;

  MatrixXd m_inv = m.inverse();
  VectorXd y = -m_inv*N.transpose()*D*x0;

  //cout << y.transpose() << endl;
  //cout << y1.transpose() << endl;

  double del_norm = (y-y1).norm();
  //cout << " dely norm = " << del_norm << endl;
  if(del_norm>1e-3){cout<<"WARNING: high error in check_contact_forces: "<<del_norm<<endl;}

  //return;
  VectorXd x = x0 + N*y;
  VectorXd x1 = x0 + N*y1;

  cout<<x.head(3).transpose()<<endl;
  cout<<x1.head(3).transpose()<<endl;

  //cout << (x-x1).transpose() << endl;
  /*cout << "y1: " << y1.transpose() << endl;
    cout << "y: " << y.transpose() << endl;*/

  cout << "cost1 = " << x1.transpose()*D*x1 << endl;
  cout << "cost = " << x.transpose()*D*x << endl; //exit(1);

}

void forcetorquesolver::solve_forcetorques(dynrecord* dynrec_, VectorXd& x, VectorXd& z){
  set_dynrec(dynrec_);
  VectorXd f; // part forces
  SpMat B; // control coupling matrix
  MatrixXd N; // null space
  solve_forcetorques_particular(x,B,f); // x is solved for zero cfs
  solve_forcetorques_null_space(B,N);

  N.conservativeResize(x.size(),N.cols());
  VectorXd y; // null space coordinates
  solve_contact_forces(x,y,N);

  MatrixXd N_cont;
  extract_N_contact(N,N_cont);
  z = -N_cont*y; // contact forces (ground reaction force)

  VectorXd delx = N*y; // change of x due to cfs
  x += delx;
  fts = x;

  return;

  cout<<fts.head(3).transpose()<<"; ";
  cout<<fts.segment(n*3,3).transpose()<<endl;
}

// constructs B matrix (for zero cfs), f vector, solves Bx=f
void forcetorquesolver::solve_forcetorques_particular(VectorXd& x, SpMat& B, VectorXd& f){
  set_B_and_f_for_zero_cfs(B,f);
  //cout<<"f:"<<endl<<f.head(10)<<endl;
  SpSolver solver;
  solver.compute(B);
  x = solver.solve(f);
}

void get_null_space(SpMat& B, int size, MatrixXd& N){
  int m = B.cols();
  SpSolver solver;
  solver.compute(B.transpose());
  if(solver.rank() < m-size){cout << "WARNING: rank shrank" << endl;}
  //cout<<"rank = "<<solver.rank()<<endl;
  
  N.resize(m,size);
  VectorXd c (m);
  c.setZero();
  for(int i=0;i<size;i++){
    if(i){c(m-i) = 0;}
    c(m-1-i) = 1;
    N.col(i) = solver.matrixQ()*c;
  }
}

// extends B to include nonzero cfs (from feet in contact)
// solves null cpace N
void forcetorquesolver::solve_forcetorques_null_space(SpMat& B, MatrixXd& N){
  int delm = 3*dynrec->get_ncontacts();
  int m = B.cols() + delm;
  B.conservativeResize(m,m); // sparse solver only worked for square matrices
  dynrec->set_forcetorque_system_contacts(B,footis);

  get_null_space(B,delm,N);

  return;
  SpMat N_sp = N.sparseView(1e-6,1e-6);
  //cout<<N_sp<<endl;
}

void submatrix_by_row_mask(MatrixXd& subM, MatrixXd& M, list<pair<int,int> >& mask, int mask_l){
  int w = M.cols();
  subM.resize(mask_l,w);
  list<pair<int,int> >::iterator it = mask.begin();
  int l = 0;
  for(;it!=mask.end();it++){
    int i = (*it).first;
    int del = (*it).second - i;
    //cout<<"i="<<i<<" del="<<del<<endl;
    subM.block(l,0,del,w) = M.block(i,0,del,w);
    l += del;
  }
}

void subvector_by_mask(VectorXd& subv, VectorXd& v, list<pair<int,int> >& mask, int mask_l){
  subv.resize(mask_l);
  list<pair<int,int> >::iterator it = mask.begin();
  int l = 0;
  for(;it!=mask.end();it++){
    int i = (*it).first;
    int del = (*it).second - i;
    subv.segment(l,del) = v.segment(i,del);
    l += del;
  }
}

void forcetorquesolver::solve_contact_forces(VectorXd& x, VectorXd& y, const MatrixXd& N){
  VectorXd c;
  set_action_penalties(c);

  DiagonalMatrix<double, Dynamic> C (c);
  MatrixXd cN = C*N;
  VectorXd cx = C*x;

  MatrixXd N0, N1;
  submatrix_by_row_mask(N0,cN,penal_mask0,mask_l0);
  submatrix_by_row_mask(N1,cN,penal_mask1,mask_l1);
  VectorXd x0, x1;
  subvector_by_mask(x0,cx,penal_mask0,mask_l0);
  subvector_by_mask(x1,cx,penal_mask1,mask_l1);

  MatrixXd N0t = N0.transpose(), N1t = N1.transpose();
  VectorXd ntx0 = N0t*x0, ntx1 = N1t*x1;
  MatrixXd ntn0 = N0t*N0, ntn1 = N1t*N1;

  double rel_error;
  int k = ntn0.cols();
  int rank0 = k;
  do {
    // zero order solution
    FullPivLU<MatrixXd> lu_decomp(ntn0);
    //lu_decomp.setThreshold(1e-13);
    while (lu_decomp.rank() > rank0) {
      lu_decomp.setThreshold(2*lu_decomp.threshold());
    }
    VectorXd y0 = lu_decomp.solve(-ntx0);
    MatrixXd Ny = lu_decomp.kernel(); // null space of ntn0
    MatrixXd Ry = lu_decomp.image(ntn0); // range of ntn0 (as well of ntn0.transpose())

    if(rank0 == lu_decomp.rank()){cout << "WARNING: decomposition threshold increased" << endl;}
    rank0 = lu_decomp.rank();
    
    VectorXd b = -(ntx1 + ntn1*y0);
    MatrixXd m (k,k);
    m << ntn1*Ny, ntn0*Ry;
    
    // first order solution
    VectorXd z = m.colPivHouseholderQr().solve(b);
    rel_error = (m*z-b).norm()/b.norm();
    rank0--;

    y = y0 + Ny*z.head(Ny.cols()); // first order lifts degeneracy of zero order solution
  } while (rel_error > 1e-6);

  return;
  check_contact_forces(x,N,c,1e-8,y,penal_mask1);
}

// cost matrix is diagonal with c squared as its diagonal
void forcetorquesolver::set_action_penalties(VectorXd& c){
  c = VectorXd::Constant(6*n,1);
  c.segment(3,3*(n-1)).setZero();

  //for(int i=3;i<c.size();i++){c(i) *= (1.1+0*sin(i));}
  for(int i=3;i<3*n;i++){c(3*n+i) = jzaxis[i];}
  if(mask_l0 == 0){cout<<"ERROR: mask0 not set"<<endl;exit(1);}
}

int mask_length(list<pair<int,int> >& mask){
  list<pair<int,int> >::iterator it = mask.begin();
  int l = 0;
  for(;it!=mask.end();it++){l += ((*it).second-(*it).first);}
  return l;
}

void forcetorquesolver::switch_torso_penalty(bool force_flag, bool torque_flag){
  penal_mask0.clear();
  vector<int> inds;
  if(force_flag){inds.push_back(0);}
  if(torque_flag){inds.push_back(3*n);}
  for(unsigned int i=0;i<inds.size();i++){
    int i0 = inds[i];
    penal_mask0.push_back(make_pair<int,int>(i0,i0+3));
  }
  mask_l0 = mask_length(penal_mask0);
  set_penal_mask1();
}

void forcetorquesolver::extract_N_contact(const MatrixXd& N, MatrixXd& N_cont){
  N_cont.resize(nf*3,N.cols());
  int k = 0;
  for(int i=0;i<nf;i++){
    for(int j=0;j<3;j++){
      N_cont.row(k++) = N.row(3*footis[i]+j);
    }
  }
}

void forcetorquesolver::set_penal_mask1(){
  penal_mask1.clear();
  list<pair<int,int> >::iterator it = penal_mask0.begin();
  int i0 = 0;
  for(;it!=penal_mask0.end();it++){
    int i1 = (*it).first;
    if(i1 > i0){
      penal_mask1.push_back(make_pair<int,int>(i0,i1));
    }
    i0 = (*it).second;
  }
  int i1 = 6*n;
  if(i1 > i0){
    penal_mask1.push_back(make_pair<int,int>(i0,i1));
  }
  mask_l1 = mask_length(penal_mask1);
}

void forcetorquesolver::set_jzaxis(){
  extvec* jzaxis_extvecs = dynrec->get_jzaxis();
  if(jzaxis == NULL){jzaxis = new double [3*n];}
  double *p = jzaxis;
  for(int i=0;i<n;i++){
    double *p1 = jzaxis_extvecs[i].get_data();
    for(int j=0;j<3;j++){*p++ = *p1++;}
  }
}

void forcetorquesolver::set_B_and_f_for_zero_cfs(SpMat& B, VectorXd& f){
  int m = 6*n;
  B.resize(m,m);
  f.resize(m); 
  dynrec->set_forcetorque_system(B,f,parentis,masses);
}

void forcetorquesolver::set_dynrec(dynrecord* dynrec_){
  dynrec = dynrec_;
  set_jzaxis();
}

void forcetorquesolver::solve_forces(dynrecord* dynrec_, VectorXd& z, VectorXd& y){
  set_dynrec(dynrec_);

  // construction of B and f
  VectorXd f; // part forces
  SpMat B; // control coupling matrix
  set_B_and_f_for_zero_cfs(B,f);

  // extention of B to include cfs from all feet
  int m = B.cols();
  int m1 = m + 3*nf;
  B.conservativeResize(m,m1);
  dynrec->set_forcetorque_system_contacts(B,footis,false);

  // extention of B and f to include torque constraints
  add_torque_constraints_to_B(B,f,z);

  // solution with zero torso force and torque  
  for(int i=0;i<2;i++){B.block(0,3*n*i,B.rows(),3) *= 0;}
  B.makeCompressed();
  SpSolver solver;
  solver.compute(B);
  VectorXd x = solver.solve(f);

  // contact forces
  y = x.tail(3*nf);
}

void forcetorquesolver::add_torque_constraints_to_B(SpMat& B, VectorXd& f, const VectorXd& z){
  int nj = z.size();
  int m = B.rows();
  B.conservativeResize(m+nj,B.cols());
  f.conservativeResize(f.size()+nj);

  list<int>::iterator it = jpart_ids.begin();
  int i = 0;
  for(;it!=jpart_ids.end();it++){
    int k = 3*(*it);
    int k1 = 3*n+k;
    for(int j=0;j<3;j++){
      B.insert(m+i,k1+j) = jzaxis[k+j];
    }
    f(m+i) = z(i);
    i++;
  }
}
