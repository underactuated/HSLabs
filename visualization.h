/////////////////////////////////////////////////
// visualization.h: declaration of functions and 
// classes used for state visualization, kinematic-
// model Open Dynamics Engine (ODE) representation
// and ODE simulations.
/////////////////////////////////////////////////
#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <vector>
#include <map>
#include <list>
#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "ode/texturepath.h"

#include "core.h"
#include "matrix.h"

void rot_ztov(dMatrix3& rot, const extvec& v);
void transpose_odematrix(dMatrix3& m);
void posrot_from_affine(dVector3& pos, dMatrix3& rot, const affine& A);
void affine_from_posrot(affine& A, const dVector3& pos, const dMatrix3& rot);
void affine_from_orientation(affine& A, const extvec* orientation);
void mod_twopi(double& a);
void euler_angles_from_affine(const affine& A, double* angles);
void evec_to_dvec(const extvec& evec, dVector3& dvec);

class kinematicmodel;
class modelplayer;
class viewpoint;
class odepart;
class trimeshmanager;

// Class visualizer provides interface to state visualization,
// ODE simulations and visualization. It stores various ODE
// constructs and ODE representation of model parts, such as
// geoms and motorized (powered) joints.
class visualizer{
  dWorldID odeworld;
  dSpaceID odespace;
  dJointGroupID contact_group;
  vector<dGeomID> geoms;
  dsFunctions fn; // drawstuff function specifying simulation loop
  const kinematicmodel* model;
  modelplayer* player;
  viewpoint* view;
  bool manual_viewpoint_flag, texture_flag;
  vector<dJointID> motors; // motorized joints
  map<dBodyID,extvec> added_forces; // for visualizing perturbing forces
  trimeshmanager* trimeshman;
  int speedup; // visualization speedup factor
public:
  visualizer();
  ~visualizer();
  inline dWorldID* get_odeworld(){return &odeworld;}
  inline dSpaceID* get_odespace(){return &odespace;}
  inline void push_geom(dGeomID geom){geoms.push_back(geom);}
  inline void set_player(modelplayer* player_){player = player_;}
  inline trimeshmanager* get_trimeshman(){return trimeshman;}
  inline const kinematicmodel* get_model() const {return model;}
  inline viewpoint* get_view(){return view;}
  inline void set_model(kinematicmodel* model_){model = model_;}
  inline int get_speedup(){return speedup;}
  void draw();
  void draw_inloop();
  void step();
  void set_viewpoint();
  void adjust_viewpoint();
  void set_flag(string flag_name, bool value);
  void simulate_odeworld(double dt_ode);
  dJointID create_contact(dContact* contact);
  void set_speedup(int f);
  void add_motor(dJointID hinge);
  void set_ode_motor_torques(const double* motor_torques);
  void get_ode_motor_angles(double* as) const;
  void get_ode_motor_adas(double* as, double* das) const;
  void get_ode_config(double* config) const;
  const odepart* get_torso_opart() const;
  void add_force(dBodyID odebody, const double* f);
private:
  void initialize_fn();
  void setup_odeworld();
  void unset_odeworld();
  void start_loop();
  void draw_forces();
};

class modelnode;

// Class odepart is a counterpart of a model node, providing access
// to modelnode's ODE representation and its state information
class odepart{
  affine A_body_geom; // relative to modelnode's body frame
  modelnode* mnode;
  dGeomID geom;
  string part_name; // corresponding body-name from xml file
  extvec capsule_to_pos; // used as foot position
  double rcap; // capsule size (for capsule geoms, 0 otherwise)
public:
  odepart(){rcap = 0;}
  inline dBodyID get_odebody() const {return dGeomGetBody(geom);}
  inline string get_part_name() const {return part_name;}
  inline modelnode* get_mnode() const {return mnode;}
  inline double get_rcap() const {return rcap;}
  void make(const xml_node<>* xnode, modelnode* mnode_, visualizer* vis);
  void get_odebody_posrot_from_body(dVector3& pos, dMatrix3& rot);
  void print(int detail_level);
  void print_ode();
  void get_com_pos(extvec& pos) const;
  void get_foot_pos(extvec& pos) const;
  void get_foot_pos(extvec& pos, bool from_body_flag) const;
  void make_fixed_joint(odepart* parent_part, visualizer* vis);
  void make_hinge_joint(odepart* parent_part, visualizer* vis);
  void get_A_ground_body_from_odebody(affine& A_ground) const;
  void get_foot_vec_ground(extvec& vec) const;
private:
  void make_ccylinder(visualizer* vis, const xml_node<>* geom_node, bool capped_flag);
  void capsule_lenposrot_from_fromto(double& len, dVector3& pos, dMatrix3& rot, const double* fromto);
  void get_A_ground_odebody(affine& A_ground) const;
};

// Class viewpoint positions camera and tracks a moving robot
// by following the torso position (using a PD controller)
class viewpoint{
  float xyz_cam_rel[3], hpr[3]; // camera params
  float xyz_ref[3], xyz_ref_rate[3]; // reference point params
  float xyz0[3], xyz_rate[3]; // torso com params
  double k0;
  bool smooth_flag; // enables PD controller
  int speedup; // visualization speedup factor
  double scale; // relative camera position scale
public:
  viewpoint();
  inline void set_smooth(bool flag_val){smooth_flag = flag_val;}
  inline void set_speedup(int speedup_){speedup = speedup_;}
  void set(const float* xyz, const float* xyz_cam, const float* hpr);
  void set(const float* xyz, const float* xyz_cam);
  void adjust(const dReal* xyz);
  void print();
  void shift_cam(float x, float y, float z);
  void set_scale(double scale_){scale = scale_;}
private:
  void hard_xyzref_update(const dReal* xyz);
  void smooth_xyzref_update(const dReal* xyz);
  void hpr_from_cam_rel();
  void set_xyz(const float* xyz, const float* xyz_cam);
};

#endif
