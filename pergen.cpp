//#include <map>
//#include "core.h"
//#include "matrix.h"
//#include "visualization.h"
//#include "model.h"
#include "lik.h"
#include "pergen.h"

// There are two ways to order limbs: 1) as in xml file (lik ordering), 
// 2) as in periodicgenerator (pergen ordering), defined below.
// Limb numbering: robot is viewed from above, moving up is moving forward.
// 0 to n/2-1: starting from top-left limb, zig-zag to the botom.
// n/2 to n-1: starting grom top-right limb, zig-zag to the botom.
// (zig-zaging means both row and column change)
// Pergen ordering is more natural for gait specification.
periodicgenerator::periodicgenerator(int n_){
  if(n_ % 2){cout << "ERROR: number of limbs must be even" << endl; exit(1);}
  n = n_;
  ts = new double [n];
  xs = new double [n];
  curvature = 0;
}

periodicgenerator::~periodicgenerator(){
  delete [] ts;
  delete [] xs;
}

// sets step duration converting 0-to-1 parameter f to t_step
// f is defined so that: f=0 <-> t_step=1/n, f=1 <-> t_step=1/2
// we use reduced time (for t_step, ts), defined as t/T
// we use reduced length (for xs), defined as x/L (?)
// xs measure foot x minus hip x
void periodicgenerator::set_step_duration(double f){
  if(f < 0 || f > 1){cout<<"ERROR: f = "<<f<<" out of bounds"<<endl;exit(1);}
  t_step = f*(1./2-1./n)+1./n;
  for(int i=0;i<2;i++){
    int jmax = n/2;
    int z = (jmax == 1)? 1 : jmax-1;
    for(int j=0;j<jmax;j++){
      int k = j+i*jmax;
      ts[k] = j*(1./2-t_step)/z + double(i)/2;
      xs[k] = ts[k] + t_step/2 - 1./2;
    }
  }
  step_duration = f;
}

// sets T, L and h (period, step_length and step_height)
void periodicgenerator::set_scales(double period_, double step_length_, double step_height_){
  period = period_;
  step_length = step_length_;
  step_height = step_height_;
}

// t here is step fraction, so it is in [0,1]
double periodicgenerator::stepx(double t){
  return (1-cos(M_PI*t))/2;
}

// see stepx comment
double periodicgenerator::stepz(double t){
  double a = sin(M_PI*t);
  return a*a;
}

// fraction of step duration
double periodicgenerator::step_frac(int limbi, double t){
  double t_lift = ts[limbi]; // liftoff moment
  if (t < t_lift) {return 0;} 
  else if (t < t_lift + t_step) {return (t-t_lift)/t_step;} 
  else {return 1;}
}

// computes foot positions in absolute coordinates, given time
void periodicgenerator::limb_positions(double time, vector<extvec>& limb_poss){
  double t = time/period; // reduced time
  int t_int = int(t);
  double t_frac = t - t_int;
  for(int i=0;i<n;i++){
    double stepf = step_frac(i,t_frac);
    double delx, delz;
    delx = (t_int + xs[i] + stepx(stepf))*step_length;
    delz = stepz(stepf)*step_height;
    extvec delpos (delx, 0, delz);
    turn_position(limb_pos0s[i],delpos,limb_poss[i]);
  }
}

void periodicgenerator::print(int detail_level){
  cout << "--- periodic generator ---" << endl;
  cout << "number of limbs = " << n << endl;
  cout << "t_step = " << t_step << endl;
  for(int i=0;i<n;i++){
    cout << "limb " << i << ":";
    cout << " ts = " << ts[i] << " xs = " << xs[i] <<endl;
  }
  cout << "period = " << period << endl;
  cout << "step_length = " << step_length << endl;
  cout << "step_height = " << step_height << endl;
  cout << "curvature = " << curvature << endl;
  if(detail_level > 0){
    cout << "limb pos0s:" << endl;
    for(int i=0;i<n;i++){limb_pos0s[i].print();}
  }
}

void periodicgenerator::change_pos0(int limbi, extvec& delpos0){
  limb_pos0s[limbi].add(delpos0);
}

void periodicgenerator::get_TLh(double TLh[3]){
  TLh[0] = period;
  TLh[1] = step_length;
  TLh[2] = step_height;
}

void periodicgenerator::set_pos0s(vector<extvec>& limb_pos0s_){
  limb_pos0s = limb_pos0s_;
  compute_max_radius();
}

void periodicgenerator::set_curvature(double curvature_){
  curvature = curvature_;
  compute_max_radius();
}

void periodicgenerator::compute_max_radius(){
  if(curvature == 0){return;}
  extvec center (0, 1./curvature, 0); // turning center
  arrayops ao (3);
  vector<extvec>::iterator it = limb_pos0s.begin();
  max_radius = 0;
  for(;it!=limb_pos0s.end();it++){
    double rad = ao.distance((*it).get_data(), center.get_data());
    if (rad > max_radius) {max_radius = rad;}
  }
}

void periodicgenerator::turn_position(extvec& pos0, extvec& delpos, extvec& pos){
  double dx, dy, dz;
  delpos.get_components(dx,dy,dz);

  if(curvature != 0){
    int s = (curvature > 0)? 1 : -1;
    double x0, y0, z0;
    pos0.get_components(x0,y0,z0);    
    double rc = 1./curvature;
    double rx = x0, ry = y0 - rc;
    double r = sqrt(rx*rx + ry*ry);
    double alpha = atan2(ry,rx);
    double beta = -s*dx/max_radius;
    double gamma = alpha - beta/2;
    double sb = 2*sin(beta/2);
    dx = r*sin(gamma)*sb;
    dy += -r*cos(gamma)*sb; 
  }
  
  extvec delpos1 (dx,dy,dz);
  delpos1.add(pos0);
  pos.copy(delpos1);
}

void periodicgenerator::get_turn_orientation(double dx, extvec* orientation){
  if(curvature != 0){
    int s = (curvature > 0)? 1 : -1;
    double psi = s*dx/max_radius;
    double rc = 1./curvature;
    orientation[0].set(rc*sin(psi), rc*(1-cos(psi)), 0);
    orientation[1].set(0, 0, psi);
  } else {
    orientation[0].set(dx, 0, 0);
    orientation[1].set(0, 0, 0);
  }
}


pergensetup::pergensetup(int n_){
  n = n_;
  pergen = new periodicgenerator(n);
  set_likpergen_map(n);
  limb_poss.resize(n);
  rec_transform.set_unity();
  rec_transform_flag = false;
}

pergensetup::~pergensetup(){
  delete pergen;
}

void pergensetup::set_TLh(double T, double L, double h){
  pergen->set_scales(T,L,h);
  v = L/T;
}

// sets record rec corresponding to time t
// record includes torso orientation (position and angles)
// and limb positions (foot positions)
// uses pergen to compute limb positions
void pergensetup::set_rec(double* rec, double t){
  extvec orientation[2];
  orientation[0].copy(torso_pos0);
  orientation[1].copy(euler_angles);
  turn_torso(t, orientation);
  orientation_to_rec(orientation, rec);  
  pergen->limb_positions(t,limb_poss);
  for(int i=0;i<n;i++){
    int j = likpergen_map[i];
    limb_poss[j].get_components(rec+6+i*3);
  }
  if(rec_transform_flag){transform_rec(rec);}
}

// maps lik indexing to pergen indexing
void pergensetup::set_likpergen_map(int n){
  switch (n) {
  case 4: { 
    int i_lik[] = {0,1,2,3};
    int i_pergen[] = {0,3,1,2};
    for(int i=0;i<n;i++){
      likpergen_map[i_lik[i]]=i_pergen[i];
    }
  } break;
  case 6: { 
    int i_lik[] = {0,1,2,3,4,5};
    int i_pergen[] = {0,3,4,1,2,5};
    for(int i=0;i<n;i++){
      likpergen_map[i_lik[i]]=i_pergen[i];
    }
  } break;
  default:
    cout << "ERROR: case " << n << " undefined in set_likpergen_map" << endl; exit(1);
  }
}

void pergensetup::set_limb_poss(int limbi, extvec& pos, double rcap){
  int j = likpergen_map[limbi];
  //pos.set_v(2,rcap+(j%2));
  pos.set_v(2,rcap);
  limb_poss[j].copy(pos);
}

void pergensetup::set_pos0s(){
  pergen->set_pos0s(limb_poss);
}

void pergensetup::set_orientation(extvec* orientation){
  torso_pos0.copy(orientation[0]);
  euler_angles.copy(orientation[1]);
}

void pergensetup::get_config_params(extvec& pos, extvec& angles, double& step_duration, double TLh[3]){
  pos.copy(torso_pos0);
  angles.copy(euler_angles);
  step_duration = pergen->get_step_duration();
  pergen->get_TLh(TLh);
}

void pergensetup::get_config_params(pgsconfigparams* pcp){
  pcp->orientation[0].copy(torso_pos0);
  pcp->orientation[1].copy(euler_angles);
  pcp->step_duration = pergen->get_step_duration();
  pergen->get_TLh(pcp->TLh);
  pcp->curvature = pergen->get_curvature();
  pcp->foot_shift = foot_shift;
}

void pergensetup::set_rec_rotation(extvec& rec_eas){
  extvec rec_transl;
  rec_transform.get_translation(rec_transl);
  set_rec_transform(rec_transl, rec_eas);
}

void pergensetup::set_rec_transform(extvec& rec_transl, extvec& rec_eas){
  extvec rec_orientation [] = {rec_transl, rec_eas};
  affine_from_orientation(rec_transform,rec_orientation);
  rec_transform_flag = true;
}

void pergensetup::transform_rec(double* rec){
  extvec orientation[2];
  rec_to_orientation(rec, orientation);
  transform_orientation(rec_transform,orientation);
  orientation_to_rec(orientation, rec);
  for(int i=0;i<n;i++){
    extvec pos0, pos;
    double* p = rec+6+i*3;
    pos0.set(p);
    rec_transform.mult(pos0,pos);
    pos.get_components(p);
  }
}

void pergensetup::copy_rec_transform(pergensetup* pgs){
  rec_transform_flag = pgs->rec_transform_flag;
  rec_transform.copy(pgs->rec_transform);
}

void pergensetup::print(int detail_level){
  cout << "---  pergen setup ---" << endl;
  cout << "number of limbs = " << n << endl;
  cout << "torso position: ";
  torso_pos0.print();
  cout << "euler angles: ";
  euler_angles.print();
  if(detail_level > 0){
    pergen->print(detail_level-1);
  }
  if(detail_level > 1){
    cout << "limb poss:" << endl;
    for(int i=0;i<n;i++){limb_poss[i].print();}
  }
}

void pergensetup::set_curvature(double curvature){
  pergen->set_curvature(curvature);
}

void pergensetup::rec_to_orientation(double* rec, extvec* orientation){
  orientation[0].set(rec);
  orientation[1].set(rec+3);
}

void pergensetup::orientation_to_rec(extvec* orientation, double* rec){
  orientation[0].get_components(rec);
  orientation[1].get_components(rec+3);
}

void pergensetup::transform_orientation(affine& A, extvec* orientation){
  affine A0, A1;
  affine_from_orientation(A0, orientation);
  A.mult(A0, A1);
  A1.get_translation(orientation[0]);
  euler_angles_from_affine(A1,orientation[1].get_data());
}

void pergensetup::turn_torso(double t, extvec* orientation){
  double tv = t*v;
  extvec torientation[2];
  pergen->get_turn_orientation(tv, torientation);
  if(torientation[1].get_v(2) == 0){
    *(orientation[0].get_data()) += tv;
  } else {
    affine A;
    affine_from_orientation(A, torientation);
    transform_orientation(A, orientation);
  }
}


pgssweeper::pgssweeper(pergensetup* pgs_, kinematicmodel* model_){
  pgs0 = pgs_;
  pgs = NULL;
  model = model_;
  pcp = new pgsconfigparams;
}

pgssweeper::~pgssweeper(){
  if(pgs){delete pgs; pgs = NULL;}
  delete pcp;
}

void pgssweeper::sweep(string param_name, double val0_, double val1, int n_val_){
  val0 = val0_;
  n_val = n_val_;
  delval = (val1-val0)/n_val;
  vali = 0;

  string sweep_names[] = {"step_duration", "period", "step_length", "step_height"};
  parami = -1;
  for(int i=0;i<4;i++){
    if(param_name == sweep_names[i]){parami = i;}
  }
  if(parami < 0){cout<<"ERROR: cannot sweep over "<<param_name<<endl;exit(1);}
  cout << "sweeping over " << param_name << ":" << endl;
}

/*bool pgssweeper::next(){return next1();
  if(vali > n_val){vali=0; return false;}
  val = val0 + vali*delval;
  vali++;

  if(pgs){delete pgs; pgs = NULL;}
  pgs0->get_config_params(pos,angles,step_duration,TLh);
  int n = pgs0->get_limb_number();
  if(parami == 0){step_duration = val;}
  else if (parami > 0 && parami < 4){TLh[parami-1] = val;}
  pgs = new pergensetup (n);
  extvec orientation[] = {pos, angles};
  setup_pergen(*pgs,orientation,step_duration);
  pgs->set_TLh(TLh);
  pgs->copy_rec_transform(pgs0);

  return true;
  }*/

bool pgssweeper::next(){
  if(vali > n_val){vali=0; return false;}
  val = val0 + vali*delval;
  vali++;

  if(pgs){delete pgs; pgs = NULL;}
  pgs0->get_config_params(pcp);
  int n = pgs0->get_limb_number();
  if(parami == 0){pcp->step_duration = val;}
  else if (parami > 0 && parami < 4){pcp->TLh[parami-1] = val;}
  pgs = new pergensetup (n);
  full_setup_pergen(*pgs, pcp);
  pgs->set_TLh(pcp->TLh);
  pgs->copy_rec_transform(pgs0);

  return true;
}

void pgssweeper::setup_pergen(pergensetup& pergensu, extvec* orientation, double step_duration_){
  model->orient_torso(orientation);
  periodicgenerator* pergen = pergensu.get_pergen();
  pergen->set_step_duration(step_duration_);
  liksolver* lik = model->get_lik();
  int n = pergensu.get_limb_number();
  assert(n == lik->get_number_of_limbs());
  double rcap = lik->get_rcap();
  setup_foot_shift(pergensu);
  for(int i=0;i<n;i++){
    extvec pos;
    lik->get_limb_pos0(i,pos);
    shift_pos0(i,pos);
    pergensu.set_limb_poss(i,pos,rcap);
  }
  pergensu.set_pos0s();
  pergensu.set_orientation(orientation);
}

void pgssweeper::full_setup_pergen(pergensetup& pergensu, pgsconfigparams* pcp0){
  pergensu.set_foot_shift(pcp0->foot_shift);
  setup_pergen(pergensu, pcp0->orientation, pcp0->step_duration);
  pergensu.set_TLh(pcp0->TLh);
  pergensu.set_curvature(pcp0->curvature);
}

void pgssweeper::setup_foot_shift(pergensetup& pergensu){
  pair<int,double> foot_shift = pergensu.get_foot_shift();
  shift_type = foot_shift.first;
  if(shift_type == 0){
    extvec shift (0, foot_shift.second, 0);
    affine* A = model->get_mnode(0)->get_A_ground();
    A->mult(shift,lat_shift);
  } else if(shift_type == 1){
    rad_shift = foot_shift.second;
  }
}

void pgssweeper::shift_pos0(int limbi, extvec& pos){
  if(shift_type < 0){return;}
  extvec delpos;
  if(shift_type == 0){
    delpos.copy(lat_shift);
    if(limbi % 2){delpos.times(-1);}
  } else if (shift_type == 1){
    double x, y, z;
    pos.get_components(x,y,z);
    double f = 1. + rad_shift/sqrt(x*x+y*y);
    pos.set(x*f,y*f,z);
  }
  pos.add(delpos);
}


pgsconfigparams::pgsconfigparams(){
  curvature = 0;
  foot_shift = pair<int,double> (-1,0);
}

void pgsconfigparams::set_TLh(double period, double step_length, double step_height){
  TLh[0] = period;
  TLh[1] = step_length;
  TLh[2] = step_height;
}
