#ifndef PLAYER_H
#define PLAYER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "matrix.h"
#include "visualization.h"
#include "model.h"

class pergensetup;
class periodic;
class configtimedertrack;
class cpccontroller;
class pgssweeper;
class heightfield;
class ghostmodel;
struct pgsconfigparams;

class modelplayer{
  kinematicmodel* model;
  int step_mode;
  bool manual_viewpoint_flag, contact_force_flag, open_loop_flag, position_control_flag, dynamics_from_simulation_flag, cpc_control_flag, record_traj_flag, torso_kicks_flag, ghost_walking_flag, clip_torque_flag;
  double play_t, play_dt;
  int config_dim, nmj; // nmj = number of motor joints
  double *play_rec, *last_motor_torques;
  pergensetup* play_pgs;
  periodic* play_per;
  configtimedertrack* tdertrack;
  double** B_transp; // transposed dynamics matrix B
  cpccontroller* play_cpc;
  vector<double*> traj_recording;
  double traj_t_limit, torque_limit;
  heightfield* hfield;
  ghostmodel* ghost;
public:
  modelplayer();
  ~modelplayer();
  inline kinematicmodel* get_model(){return model;}
  void load_model(string fname);
  void step();
  void test(int testi);
  void set_flag(string flag_name, bool value);
  void set_jangles_with_lik(double* rec);
  void print_limb_pos0s();
  void orient_torso(extvec* orientation);
  void setup_pergen(pergensetup& pergensu, extvec* orientation, double step_duration);
  pergensetup* make_pergensu(string config_fname, int setup_id);
  void full_setup_pergen(pergensetup& pergensu, pgsconfigparams pcp);
  //void get_pgs_config_params(string& rec_str, string& fname, extvec* orientation, double& step_duration, double& period, double& step_length, double& step_height);
  void get_pgs_config_params(string& rec_str, pgsconfigparams& confparams);
  void play_pergensu(pergensetup* pgs);
  void set_play_dt(double dt){play_dt = dt;}
  void prepare_per_traj_dyn(periodic& per, pergensetup* pgs, int n_t);
  double measure_cot(pergensetup* pgs, int n_t);
  void pergensu_config_string(pergensetup* pgs, string& str);
  void measure_cot_sweep(pergensetup* pgs, int n_t, string param_name, double val0, double val1, int n_val);
  void test_dynamics(pergensetup* pgs); // probably temporary
  void simulate_pergensu(pergensetup* pgs, double t0);
  void speedup_draw(int f);
  void open_loop_test(pergensetup* pgs, double t0);
  void position_control_test(pergensetup* pgs, double t0);
  void cpc_test(pergensetup* pgs, double t0);
  void record_pos_control_traj(pergensetup* pgs, double t0, double traj_duration);
  void load_draw_model(string fname); // prob tmp
  void record_per_traj(pergensetup* pgs);
  void record_per_traj_sweep(pergensetup* pgs, string param_name, double val0, double val1, int n_val);
  void uneven_ground_test();
  void set_torque_limit(double torque);
  void test_lik_solvers();
  void shift_view(double x, double y, double z);
private:
  inline visualizer* get_vis(){return model->get_vis();}
  void test0();
  void test1();
  void test2();
  void test3();
  void test4();
  void get_rec_str(string& rec_str, string fname, int rec_id);
  void play_pergensu();
  void make_play_rec();
  void simulate_ode();
  void set_open_loop_torques();
  void init_play_config(pergensetup* pgs);
  void setup_per_controller(pergensetup* pgs, double t0);
  void unset_per_controller();
  void set_position_control_torques();
  void estimate_B();
  void save_last_motor_torques(double* torques);
  void setup_cpc_controller();
  void set_cpc_torques();
  void add_traj_record();
  void save_traj_recording(string fname, int rec_len);
  void estimate_B_by_regression();
  void kick_torso(); // experimental
  void set_default_flags();
  void torso_velocity(); // experimental
  dBodyID get_torso_odebody();
  void fall_check(double hc);
  void linear_feedback_control(double* torques, double** x0, double** x, double k1, double k2);
  void set_ghost_cpc_state(double** ders);
  void max_motor_torque(); // experimental
  void set_ode_motor_torques(double* motor_torques);
  void clip_torque(double* motor_torques);
  void check_model_loaded();
};

#endif
