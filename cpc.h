#ifndef CPC_H
#define CPC_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "matrix.h"
#include "model.h"
#include "periodic.h"

using namespace Eigen;

struct cand{
  int tpi;
  double t0, s, loss;
  cand(int tpi, double t0, double s, double loss);
  cand(){};
  void print();
};

class efficientdata;

class cpccontroller{
  int q_dim, chi_dim, psi_dim; // q_dim = config_dim, chi_dim = nmj = actuated dim, psi_dim = q_dim - chi_dim = unactuated dim
  MatrixXd B, Bt, B_chi, b, bbt, btil0; // Bt = B transposed
  double** target_points;
  int tps_size; // tps = target point set
  double* current_state;
  double w, sg, k0, kc, tauc;
  int n_cand;
  ColPivHouseholderQR<MatrixXd> B_chi_dec;
  cand last_cand;
  bool poscontrol_flag, goal_t0s_flag, tpdist_switch_flag;
  int poscontrol_tpi;
  efficientdata* effdata;
  vector<int> mask;
  bool low_tpdist;
public:
  cpccontroller(const kinematicmodel* model);
  ~cpccontroller();
  void set_target_points_by_per(periodic* per);
  void set_B_by_player(double** B_transp); // maybe call it set_B_transp ?
  void print();
  void print(int detail_level);
  void print_target_points();
  void set_current_state(const double* config, const double* dconfig);
  void set_w_sg_nd_k0_kc_tauc(double w, double sg, int nd, double k0, double kc_, double tauc_);
  void get_motor_torques(double* torques);
  void load_tpset(string fname);
void get_torques0(double* torques); // temporaty/experimental
  void set_flag(string flag_name, bool value);
  void setup_effdata();
private:
  void compute_b();
  void compute_prox_loss_over_tpset(); // probably temporary
  void state_to_qdq(double* state, VectorXd& q, VectorXd& dq);
  void state_to_chidchi(double* state, VectorXd& chi, VectorXd& dchi);
  void compute_prox_loss_over_tps(const list<int>& tpis, map<double,cand>& loss_cands);
  double prox_loss(double t0, double s, double cc0, double cc1);
  void candidates(list<cand>& cands);
  void controls(double* torques, double k, const cand& can);
  void controls(VectorXd& tau, double k, const cand& can);
  void controls1(VectorXd& tau, VectorXd& del_tau, double k, const cand& can);
  void get_tau(int tpi, VectorXd& tau);
  int cost(double k, const list<cand>& cands, MatrixXd& cand_taus);
  int cost2(double k, const list<cand>& cands, MatrixXd& cand_taus);
  int cost1(double k, const list<cand>& cands, MatrixXd& cand_taus);
  void get_poscontrol_torques(double* torques);
  void tpis_subset(list<int>& tpis);
  void tpis_effdata(list<int>& tpis, int n);
  void apply_mask(double* state);
  void tp_dist_check(); // experim temp
};

void fit_plane(extvec& plane, list<extvec>& points);

#endif
