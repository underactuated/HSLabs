/////////////////////////////////////////////////
// ghost.h: declaration of ghostmodel class
// implementing the ghost walking method enabling
// walking on various types of strongly uneven
// terrains. This is a very low bandwidth method
// in terms of computational and sensory demands.
/////////////////////////////////////////////////
#ifndef GHOST_H
#define GHOST_H

#include "matrix.h"
#include "visualization.h"
#include "model.h"
#include "lik.h"
#include "odestate.h"
#include "geom.h"

// Ghost walking method summary: In ghost walking
// approach the current state is mapped on to 
// a ghost state, corresponding to walking on flat
// surface. All the torques are then computed with
// ordinary position control using a single walking
// cycle (on flat terrain) as the target trajectory
// for the ghost. No optimization or planning is
// involved in this approach. The only information
// used by the mapping function is the direction of
// gravity and (instantaneous) clearance under
// every foot.
/////////////////////////////////////////////////
// Class ghostmodel impements the ghost model
// counterpart of the ghost model approach and
// provides methods for mapping between the original
// model and its ghost. 
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
