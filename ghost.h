#ifndef GHOST_H
#define GHOST_H

#include "matrix.h"
#include "visualization.h"
#include "model.h"
#include "lik.h"
#include "odestate.h"
#include "geom.h"

class ghostmodel{
  const visualizer* vis;
  const heightfield* hfield;
  int config_dim, nmj, nfeet;
  double dt, rcap, foot_min_z;
  double *config, *gconfig, *lik_rec; // (g)config - (g)model config
  const odepart* torso_opart;
  vector<odepart*> foot_oparts;
  kinematicmodel* gmodel; // gmodel = ghost model
  vector<extvec> limb_poss;
  configtimedertrack* tdertrack;
  int idle;
  extvec torso_pos, torso_eas, torso_com, parent_pos; // maybe change func names to orient
  extvec surf_normal, plane_zf;
  affine surf_rot;
  bool horizontal_flag, adaptive_orientation_flag;
  double** acd; // array of config_dim-sized elements
  const liksolver *lik, *glik; // (g)model LIK
public:
  ghostmodel(const visualizer* vis, const heightfield* hfield, double dt);
  ~ghostmodel();
  kinematicmodel* nonvis_clone(const kinematicmodel* model);
  void get_motor_adas(double* as, double* das);
  void set_torso_feet_oparts(const kinematicmodel* model);
  void set_surf_rot(const extvec& eas);
  void rotate_surf();
private:
  void shift_limb_poss();
  void set_lik_rec();
  void set_plane_zf();
  void surf_rotate_pos(extvec& pos);
  void surf_rotate_torso_pos();
  void set_torso_pos();
  void set_surf_rot_from_normal();
  void set_gmodel_limb_bends();
};

#endif
