/////////////////////////////////////////////////
// periodic.h: declaration of the class periodic,
// dealing with periodic trajectories. Presented
// with a sequence of configuration points (and
// time step between the points), representing a
// periodic trajectory, it fills in the missing
// information, such as velocities and accelerations
// of all parts, as well as forces and torques.
// It should eventually be used for computing
// objective function derivatives.
/////////////////////////////////////////////////
#ifndef PERIODIC_H
#define PERIODIC_H

#include "core.h"
#include "model.h"
#include "pergen.h"

class dynpart;
class dynrecord;
class forcetorquesolver;

// Class periodic contains information determining
// model dynamics. Periodic computes time derivatives
// of quantities of interest along the trajectory,
// as well as actuations that realize it.
class periodic{
  const kinematicmodel* model;
  vector<dynpart*> dynparts; // model part info determining dynamics
  int n_t, config_dim, traj_size, nmj; // traj_size is slightly larger than n_t to simplify computation of derivatives over one period
  double** traj; // config traj, no vels or torques
  double dt_traj, rcap; // time step and foot cap radius
  dynrecord** dynrecs; // dynamic quantities along the trajectory
  int *parentis, *footis; // parent ids (enumerate by child ids), foot node ids (enumerated by limb ids in lik-order)
  double* masses; // part masses
  int nfeet; // number of feet
  double **vel_traj, **computed_torques; // velocities and torques along traj
  forcetorquesolver* ftsolver; // force and torque solver
  double min_cfz, max_mu; // minimum ground reaction force, maximum friction coeff
public:
  periodic(const kinematicmodel* model);
  ~periodic();
  inline dynpart* get_dynpart(int i) const {return dynparts[i];}
  inline double* get_masses() const {return masses;}
  inline int get_number_of_dynparts() const {return dynparts.size();}
  inline int* get_parentis() const {return parentis;}
  inline int* get_footis() const {return footis;}
  inline int get_nfeet() const {return nfeet;}
  inline int get_nt(){return n_t;}
  inline void get_contforce_stat(double* stat){stat[0] = min_cfz; stat[1] = max_mu;}
  inline double* get_computed_torques(int i) const {return computed_torques[i % n_t];}
  void print();
  void record_trajectory(pergensetup* pgs, int n_t);
  void print_trajectory();
  void print_trajectory(int i);
  void print_trajectory(int i, int m);
  void compute_dynrecs();
  void print_dynrecs(int id);
  void print_dynrecs(int id, int i, int m);
  void print_dynrec(int i);
  void recompute_dynparts();
  void compute_dynrec_ders();
  void switch_torso_penalty(bool force, bool torque);
  void check_solve_ft(); // temporary
  void hinge_joint_part_ids(list<int>& ids) const;
  void compute_vel_traj();
  double work_over_period();
  double get_total_mass();
  void get_motor_torques(double* motor_torques);
  void solve_torques_contforces(int i, double* torques, double* contforces);
  void solve_contforces_given_torques(int i, double* contforces, double* torques);
  void compute_torques_over_period();
  void get_motor_adas(int tsi, double* as, double* das);
  void get_complete_traj_rec(int tsi, double* rec);
  void get_complete_traj(double** complete_traj);
private:
  void set_dynparts();
  void clear_traj();
  void new_dynrecs();
  void delete_dynrecs();
  void set_footset(set<modelnode*>& foot_set);
  void new_vel_traj();
  void delete_vel_traj();
  void analyze_contforces(const double* contforces);
  void new_computed_torques();
  void delete_computed_torques();
};

#endif
