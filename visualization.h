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

void rot_ztov(dMatrix3& rot, extvec& v);
void transpose_odematrix(dMatrix3& m);
void posrot_from_affine(dVector3& pos, dMatrix3& rot, affine& A);
void affine_from_posrot(affine& A, const dVector3& pos, const dMatrix3& rot);
void affine_from_orientation(affine& A, extvec* orientation);
void mod_twopi(double& a);
void euler_angles_from_affine(affine& A, double* angles);

// ##### VISUALIZER #####
// includes minimal functionality from ode and ds to visualize robot state

class kinematicmodel;
class modelplayer;
class viewpoint;
class odepart;
class trimeshmanager;

class visualizer{
  dWorldID odeworld;
  dSpaceID odespace;
  dJointGroupID contact_group;
  vector<dGeomID> geoms;
  dsFunctions fn;
  kinematicmodel* model;
  modelplayer* player;
  viewpoint* view;
  bool manual_viewpoint_flag, texture_flag;
  vector<dJointID> motors;
  map<dBodyID,extvec> added_forces;
  trimeshmanager* trimeshman;
public:
  visualizer();
  ~visualizer();
  inline dWorldID* get_odeworld(){return &odeworld;}
  inline dSpaceID* get_odespace(){return &odespace;}
  inline void push_geom(dGeomID geom){geoms.push_back(geom);}
  inline void set_player(modelplayer* player_){player = player_;}
  inline trimeshmanager* get_trimeshman(){return trimeshman;}
  inline kinematicmodel* get_model(){return model;}
  inline viewpoint* get_view(){return view;}
  inline void set_model(kinematicmodel* model_){model = model_;}
  void initialize_fn();
  void setup_odeworld();
  void unset_odeworld();
  void start_loop();
  void draw();
  void draw_inloop();
  void step();
  void set_viewpoint();
  void adjust_viewpoint();
  void set_flag(const string flag_name, const bool value);
  void simulate_odeworld(const double dt_ode);
  dJointID create_contact(dContact* contact);
  void set_speedup(const int f);
  void add_motor(dJointID hinge);
  void set_ode_motor_torques(const double* motor_torques);
  void get_ode_motor_angles(double* as);
  void get_ode_motor_adas(double* as, double* das);
  void get_ode_config(double* config);
  odepart* get_torso_opart();
  void add_force(dBodyID body, const double* f);
  void draw_forces();
private:
  void trimesh_test(); // temp
};

class modelnode;

class odepart{
  affine A_geom; // in modelnode's body frame
  modelnode* mnode;
  dGeomID geom;
  string part_name;
  extvec capsule_to_pos; // used as foot position
  double rcap; // capsule size (for capsule geoms, 0 otherwise)
public:
  odepart(){rcap = 0;}
  inline dBodyID get_body(){return dGeomGetBody(geom);}
  inline string get_part_name(){return part_name;}
  inline modelnode* get_mnode(){return mnode;}
  inline double get_rcap(){return rcap;}
  void make(const xml_node<>* xnode, modelnode* mnode_, visualizer* vis);
  void capsule_lenposrot_from_fromto(double& len, dVector3& pos, dMatrix3& rot, const double* fromto);
  void get_body_posrot_from_frame(dVector3& pos, dMatrix3& rot);
  void print(const int detail_level);
  void print_ode();
  void get_com_pos(extvec& pos);
  void get_foot_pos(extvec& pos);
  void get_foot_pos(extvec& pos, const bool from_body_flag);
  void make_fixed_joint(odepart* parent_part, visualizer* vis);
  void make_hinge_joint(odepart* parent_part, visualizer* vis);
  void get_frame_A_ground_from_body(affine& A_ground);
private:
  void get_ode_body_A_ground(affine& A_ground);
  void make_ccylinder(visualizer* vis, const xml_node<>* geom_node, const bool capped_flag);
};

class viewpoint{
  float xyz_ref[3], xyz_cam_rel[3], hpr[3];
  float k0, xyz0[3], xyz_rate[3], xyz_ref_rate[3];
  bool smooth_flag;
public:
  viewpoint();
  inline void set_smooth(const bool flag_val){smooth_flag = flag_val;}
  void set(const float* xyz, const float* xyz_cam, const float* hpr);
  void set(const float* xyz, const float* xyz_cam);
  void adjust(const dReal* xyz);
  void print();
  void shift_cam(const float x, const float y, const float z);
private:
  void hard_xyzref_update(const dReal* xyz);
  void smooth_xyzref_update(const dReal* xyz);
  void hpr_from_cam_rel();
  void set_xyz(const float* xyz, const float* xyz_cam);
};

#endif
