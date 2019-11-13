/////////////////////////////////////////////////
// pergen.h: declaration of classes tasked with
// generating periodic gates. There are two main
// classes: periodicgenerator (pergen for short)
// and pergen setup. periodicgenerator contains
// minimal information needed to generate periodic
// motion of feet, possible along a curved path.
// (pergen has no information about torso).
// pergensetup contains additional information
// needed to generate full robot gate (including
// its torso).
/////////////////////////////////////////////////
#ifndef PERGEN_H
#define PERGEN_H

#include "matrix.h"
#include "model.h"

// Class periodicgenerator generates peiodic motion
// of feet. Periodic pattern is formulated in
// reduced time and distance, measured in period
// and step length respectivly. When curvature = 0,
// walking motion is along x axis. Walking while
// turning is produced by transforming x-y plane so
// that strait lines of constant y transform into
// concentric circles.
class periodicgenerator{
  int n; // number of limbs
  double t_step; // step duration in reduced time
  double *ts, *xs; // reduced time and space liftoff moments
  double period, step_length, step_height;
  vector<extvec> limb_pos0s; // default foot positions (pergen ordering)
  double step_duration; // parameter in [0,1]
  double curvature, max_radius; // path curvature, max turning-center to hip distance
public:
  periodicgenerator(int n);
  ~periodicgenerator();
  inline double get_period() const {return period;}
  inline double get_step_length(){return step_length;}
  inline double get_step_duration(){return step_duration;}
  inline double get_curvature(){return curvature;}
  void set_step_duration(double f); // f is in [0,1]
  void set_scales(double period, double step_length, double step_height);
  void limb_positions(double time, vector<extvec>& limb_poss);
  void set_pos0s(const vector<extvec>& limb_pos0s);
  void change_pos0(int limbi, const extvec& delpos0);
  void print(){print(0);}
  void print(int detail_level);
  void get_TLh(double TLh[3]);
  void set_curvature(double curvature);
  void get_turn_orientation(double dx, extvec* orientation);
private:
  double stepx(double t) const;
  double stepz(double t) const;
  double step_frac(int limbi, double t) const;
  void compute_max_radius();
  void turn_position(const extvec& pos0, const extvec& delpos, extvec& pos);
};

struct pgsconfigparams;

// Class pergensetup ( = periodic generator setup) completes
// periodic cycle specification by incrorporating torso
// information. It can also additionally transform
// trajectory (record) to generate motion in arbitrary
// direction. (recall that pergen only generate motion
// in x direction)
class pergensetup{
  int n; // number of limbs
  periodicgenerator* pergen;
  map<int,int> likpergen_map; // LIK-pergen correspondence 
  vector<extvec> limb_poss; // limb positions = foot positions (pergen ordering)
  double v; // velocity = step_length/period
  extvec torso_pos0, euler_angles; // torso orientation
  affine rec_transform; // trajectory transformation
  bool rec_transform_flag;
  pair<int,double> foot_shift; // foot shift parameter = pair<type,value>
public:
  pergensetup(int n);
  ~pergensetup();
  inline periodicgenerator* get_pergen() const {return pergen;}
  inline int get_limb_number() const {return n;}
  inline int get_config_dim() const {return 6+3*n;}
  inline double get_period() const {return get_pergen()->get_period();}
  inline pair<int,double> get_foot_shift(){return foot_shift;}
  inline void set_TLh(const double TLh[3]) {set_TLh(TLh[0],TLh[1],TLh[2]);}
  inline void set_foot_shift(const pair<int,double> foot_shift_){foot_shift = foot_shift_;}
  void set_TLh(double T, double L, double h);
  void set_rec(double* rec, double t);
  void set_limb_poss(int limbi, const extvec& pos, double rcap);
  void set_pos0s();
  void set_orientation(const extvec* orientation);
  void get_config_params(extvec& pos, extvec& angles, double& step_duration, double TLh[3]);
  void get_config_params(pgsconfigparams* pcp) const;
  void set_rec_rotation(const extvec& rec_eas);
  void copy_rec_transform(const pergensetup* pgs);
  void print(){print(0);}
  void print(int detail_level);
  void set_curvature(double curvature);
private:
  void set_likpergen_map(int n);
  void set_rec_transform(const extvec& rec_transl, const extvec& rec_eas);
  void transform_rec(double* rec);
  void rec_to_orientation(const double* rec, extvec* orientation);
  void orientation_to_rec(const extvec* orientation, double* rec);
  void transform_orientation(const affine& A, extvec* orientation);
  void turn_torso(double t, extvec* orientation);
};

// Class pgssweeper iterates over pgs configurations by varying
// one of the config parameters.
class pgssweeper{
  const pergensetup *pgs0; // reference pgs
  pergensetup *pgs;
  int parami, n_val, vali; // parameter index, number of values, value index
  double val0, delval, val; // initial value, increment, current value
  const kinematicmodel* model;
  pgsconfigparams* pcp;
  int shift_type;
  extvec lat_shift; // lateral shift value (also direction)
  double rad_shift; // radial shift value
public:
  pgssweeper(const pergensetup* pgs, const kinematicmodel* model);
  ~pgssweeper();
  inline pergensetup* get_pgs() const {return pgs;}
  inline double get_val(){return val;}
  void sweep(string param_name, double val0, double val1, int n_val);
  bool next();
  void partial_setup_pergen(pergensetup& pergensu, const extvec* orientation, double step_duration);
  void setup_pergen(pergensetup& pergensu, const pgsconfigparams* pcp);
private:
  void shift_pos0(int limbi, extvec& pos);
  void setup_foot_shift(pergensetup& pergensu);
};

// Structure pgsconfigparams stores pergen-setup config parameters.
struct pgsconfigparams{
  string fname;
  extvec orientation[2];
  double step_duration;
  double TLh[3];
  double curvature;
  pair<int,double> foot_shift; // 0 lateral shift, 1 radial shift
  pgsconfigparams();
  void set_TLh(double period, double step_length, double step_height);
};

#endif
