#include <vector>
#include <map>
#include <list>
#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include "ode/texturepath.h"

void rot_ztov(dMatrix3& rot, extvec& v);
void transpose_odematrix(dMatrix3& m);
void posrot_from_affine(dVector3& pos, dMatrix3& rot, affine& A);
void affine_from_posrot(affine& A, const dVector3& pos, const dMatrix3& rot);
void affine_from_orientation(affine& A, extvec* orientation);
void mod_twopi(double& a);
void euler_angles_from_affine(affine& A, double* angles);

// ##### VISUALIZER #####
// includes minimal functionality from ode and ds to visualize robot state

struct kinematicmodel;
struct modelplayer;
struct viewpoint;
struct odepart;
struct trimeshmanager;

struct visualizer{
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
  visualizer();
  ~visualizer();
  void initialize_fn();
  void setup_odeworld();
  void unset_odeworld();
  void start_loop();
  void draw();
  void draw_inloop();
  dWorldID* get_odeworld(){return &odeworld;}
  dSpaceID* get_odespace(){return &odespace;}
  void push_geom(dGeomID geom){geoms.push_back(geom);}
  void step();
  void set_player(modelplayer* player_){player = player_;}
  void set_viewpoint();
  void adjust_viewpoint();
  void set_flag(string flag_name, bool value);
  void simulate_odeworld(double dt_ode);
  dJointID create_contact(dContact* contact);
  void set_speedup(int f);
  void add_motor(dJointID hinge);
  void set_ode_motor_torques(double* motor_torques);
  void get_ode_motor_angles(double* as);
  void get_ode_motor_adas(double* as, double* das);
  void get_ode_config(double* config);
  odepart* get_torso_opart();
  void add_force(dBodyID body, double* f);
  void draw_forces();
  trimeshmanager* get_trimeshman(){return trimeshman;}
  kinematicmodel* get_model(){return model;}
  viewpoint* get_view(){return view;}
private:
  void trimesh_test(); // temp
};

struct modelnode;

struct odepart{
  affine A_geom; // in modelnode's body frame
  modelnode* mnode;
  dGeomID geom;
  string part_name;
  extvec capsule_to_pos; // used as foot position
  void make(xml_node<>* xnode, modelnode* mnode_, visualizer* vis);
  void capsule_lenposrot_from_fromto(double& len, dVector3& pos, dMatrix3& rot, double* fromto);
  dBodyID get_body(){return dGeomGetBody(geom);}
  void get_body_posrot_from_frame(dVector3& pos, dMatrix3& rot);
  void print(int detail_level);
  void print_ode();
  string get_part_name(){return part_name;}
  modelnode* get_mnode(){return mnode;}
  void get_com_pos(extvec& pos);
  void get_foot_pos(extvec& pos);
  void get_foot_pos(extvec& pos, bool from_body_flag);
  void make_fixed_joint(odepart* parent_part, visualizer* vis);
  void make_hinge_joint(odepart* parent_part, visualizer* vis);
  void get_frame_A_ground_from_body(affine& A_ground);
private:
  void get_ode_body_A_ground(affine& A_ground);
  void make_ccylinder(visualizer* vis, xml_node<>* geom_node, bool capped_flag);
};

struct viewpoint{
  float xyz_ref[3], xyz_cam_rel[3], hpr[3];
  float k0, xyz0[3], xyz_rate[3], xyz_ref_rate[3];
  bool smooth_flag;
  viewpoint();
  void set(float* xyz, float* xyz_cam, float* hpr);
  void set(float* xyz, float* xyz_cam);
  void adjust(const dReal* xyz);
  void print();
  void set_smooth(bool flag_val){smooth_flag = flag_val;}
  void shift_cam(float x, float y, float z);
private:
  void hard_xyzref_update(const dReal* xyz);
  void smooth_xyzref_update(const dReal* xyz);
  void hpr_from_cam_rel();
  void set_xyz(float* xyz, float* xyz_cam);
};

