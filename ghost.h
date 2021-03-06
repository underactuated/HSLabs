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
  const visualizer* vis; // visualizer (for accessing original model)
  const heightfield* hfield; // uneven terrain
  int config_dim, nmj, nfeet; // configuration dimension, number of motorized joints, number of feet
  double dt, rcap, foot_min_z; // time step, foot cap radius, minimum foot height
  double *config, *gconfig, *lik_rec; // (g)config = (g)model config
  const odepart* torso_opart; // model (= original model) torso odepart
  vector<odepart*> foot_oparts; // model foot odeparts
  kinematicmodel* gmodel; // gmodel = ghost model
  vector<extvec> limb_poss;
  configtimedertrack* tdertrack;
  int idle;
  extvec torso_pos, torso_eas, torso_com, parent_pos; // torso orientaion (_pos + _eas); torso com position in ground frame; torso joint A_parent's translation 
  extvec surf_normal; // surface normal = fitted plane normal
  extvec fitted_plane, torso_plane; // fitted plane, torso plane (parallel to fitted plane, contains COM)
  affine surf_rot; // rotation around COM that takes surf_normal to (0,0,1)
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
  void transform_limb_poss();
  void set_lik_rec();
  void set_torso_plane();
  void surf_rotate_pos(extvec& pos);
  void surf_rotate_torso_pos();
  void set_torso_pos();
  void set_surf_rot_from_normal();
  void set_gmodel_limb_bends();
  void compute_surf_normal(list<extvec>& points);
  double get_plane_z(extvec& plane, extvec& pos);
};

#endif
