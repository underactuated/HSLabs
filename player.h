/////////////////////////////////////////////////
// player.h: declaration of modelplayer class. It
// serves as a main testbed for this project. Player
// contains numerous tests, can play model gaits
// and run simulations. It is actively evolving.
// For now modelplayer is intended to provide necessary 
// interface for all the available functionality.
// main.cpp contains many examples.
/////////////////////////////////////////////////
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

// Class modelplayer has a number of flags for configuring the player:
// manual_viewpoint_flag - viewpoint is adjusted manually with a mouse
// 	in the ds visualization window. 
// contact_force_flag - to print contact force statistics during COT
// 	measurements.
// open_loop_flag - open loop (torque) controller
// position_control_flag - (open loop) position controller
// dynamics_from_simulation_flag - to estimate dynamics matrix B from ODE
//	simulation. B is needed for CPC. This is tempropy, eventually
//	B estimation should be model based.
// cpc_control_flag - configuration path control (CPC) controller
// record_traj_flag - to record simulated trajectory.
// torso_kicks_flag - to subject torso to "kicks".
// ghost_walking_flag - to use ghost walking method.
// clip_torque_flag - to limit motor torques.
class modelplayer{
  kinematicmodel* model;
  int step_mode; // mode of modelplayer::step function, called from visualizer loop
  bool manual_viewpoint_flag, contact_force_flag, open_loop_flag, position_control_flag, dynamics_from_simulation_flag, cpc_control_flag, record_traj_flag, torso_kicks_flag, ghost_walking_flag, clip_torque_flag;
  double play_t, play_dt; // playing time and time step size
  int config_dim, nmj; // nmj = number of motor joints
  double *play_rec, *last_motor_torques; // play_rec - for set_jangles_with_lik(rec) 
  pergensetup* play_pgs; // periodic generator setup
  periodic* play_per; // object for computing and managing periodic trajectories
  configtimedertrack* tdertrack; // for tracking configuration time derivatives
  double** B_transp; // transposed dynamics matrix B
  cpccontroller* play_cpc;
  vector<double*> traj_recording; // storage for recorded simulated trajectory
  double traj_t_limit, torque_limit; // recording limit, torque limit
  heightfield* hfield; // for uneven terrain
  ghostmodel* ghost; // for ghost walking
public:
  modelplayer();
  ~modelplayer();
  inline kinematicmodel* get_model() const {return model;}
  inline void set_play_dt(double dt){play_dt = dt;}
  void load_model(string fname);
  void step();
  void test(int testi);
  void set_flag(string flag_name, bool value);
  void set_jangles_with_lik(const double* rec);
  void print_hip_poss();
  void orient_torso(const extvec* orientation) const;
  void partial_setup_pergen(pergensetup& pergensu, const extvec* orientation, double step_duration);
  pergensetup* make_pergensu(string config_fname, int setup_id);
  void setup_pergen(pergensetup& pergensu, const pgsconfigparams& pcp);
  void get_pgs_config_params(const string& rec_str, pgsconfigparams& confparams);
  void play_pergensu(pergensetup* pgs);
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
  inline visualizer* get_vis() const {return model->get_vis();}
  void test0();
  void test1();
  void test2();
  void test3();
  void test4();
  void get_rec_str(string& rec_str, string fname, int rec_id);
  void play_pergensu_step();
  void make_play_rec();
  void simulate_ode();
  void set_open_loop_torques();
  void init_play_config(pergensetup* pgs);
  void setup_per_controller(pergensetup* pgs, double t0);
  void unset_per_controller();
  void set_position_control_torques();
  void estimate_B();
  void save_last_motor_torques(const double* torques);
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
