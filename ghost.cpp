//#include "core.h"
//#include "matrix.h"
//#include "visualization.h"
//#include "model.h"
//#include "lik.h"
//#include "odestate.h"
//#include "geom.h"
#include "ghost.h"

void fit_plane(extvec& plane, list<extvec>& points);

ghostmodel::ghostmodel(visualizer* vis_, heightfield* hfield_, double dt_){
  vis = vis_;
  hfield = hfield_;
  dt = dt_;
  kinematicmodel* model = vis->get_model();
  config_dim = model->get_config_dim();
  nmj = model->number_of_motor_joints();
  //ao.set_n(config_dim);
  acd = new_2d_array(3,config_dim);
  config = acd[0];
  gconfig = acd[1];
  lik_rec = acd[2];
  set_torso_feet_oparts(model);
  gmodel = nonvis_clone(model);
  tdertrack = new configtimedertrack (2,config_dim,dt);
  idle = 3;
  surf_normal.set(0,0,1);
  surf_rot.set_unity();
  rotate_surf();
  set_plane_zf();
  horizontal_flag = true;
  adaptive_orientation_flag = true;//false;
  lik = model->get_lik();
  glik = gmodel->get_lik();
  rcap = lik->get_rcap();
}

ghostmodel::~ghostmodel(){
  delete_2d_array(acd,3);
  delete gmodel;
  delete tdertrack;
}

// clones a model without ode representation in visualizer
kinematicmodel* ghostmodel::nonvis_clone(const kinematicmodel* model){
  kinematicmodel* clone = new kinematicmodel (false);
  clone->load_fromxml(model->get_xmlfname());
  clone->recompute_modelnodes();
  return clone;
}

// motor adas of the ghost model
void ghostmodel::get_motor_adas(double* as, double* das){
  rotate_surf();
  shift_limb_poss();
  set_lik_rec(); //print_array<double>(lik_rec,config_dim);cout<<endl;
  gmodel->set_jvalues_with_lik(lik_rec);
  gmodel->get_jvalues(gconfig);
  tdertrack->push_config(gconfig);
  double** ders = tdertrack->get_ders();
  arrayops ao (nmj);
  if(idle){idle--; return;}
  ao.assign(as,ders[0]+6);
  ao.assign(das,ders[1]+6);
}

// gets torso and feet odeparts from model
void ghostmodel::set_torso_feet_oparts(kinematicmodel* model){
  torso_opart = model->get_vis()->get_torso_opart();
  torso_opart->get_mnode()->get_joint()->get_A_parent()->get_translation(parent_pos);

  set<modelnode*> foot_mnodes;
  model->get_foot_mnodes(foot_mnodes);
  vector<odepart*>* oparts = model->get_odeparts();
  vector<odepart*>::iterator it = oparts->begin();
  for(;it!=oparts->end();it++){
    modelnode* mnode = (*it)->get_mnode();
    if(foot_mnodes.count(mnode)){
      foot_oparts.push_back(*it);
    }
  }
  nfeet = foot_oparts.size();
}

// transforms model foot positions into ghost foor positions
void ghostmodel::shift_limb_poss(){
  double zmin = 1e10;
  list<extvec> points;
  limb_poss.resize(nfeet);
  for(int i=0;i<nfeet;i++){
    extvec* lp = &(limb_poss[i]);
    foot_oparts[i]->get_foot_pos(*lp, true);
    double x, y, z;
    lp->get_components(x,y,z);
    surf_rotate_pos(*lp);
    double z1 = lp->get_v(2);
    if(z1 < zmin){zmin = z1;}
    double h = hfield->get_h(x,y);
    h += rcap*(1./surf_normal.get_v(2)-1.);
    lp->set_v(2,z-h);
    extvec point (x,y,h+1.);
    point.subtract(torso_com);
    points.push_back(point);
  }
  foot_min_z = zmin;
  if(adaptive_orientation_flag){
    fit_plane(surf_normal,points);
  }
  surf_normal.normalize();
  set_surf_rot_from_normal();
  //cout<<acos(surf_normal.get_v(2))*180/M_PI<<endl;
}

// sets lik_rec for the ghost
void ghostmodel::set_lik_rec(){
  vis->get_ode_config(config);
  set_gmodel_limb_bends();
  surf_rotate_torso_pos();
  torso_pos.get_components(lik_rec); //torso_pos.print();
  double* p = lik_rec+6;
  vector<extvec>::iterator it = limb_poss.begin();
  for(;it!=limb_poss.end();it++){
    (*it).get_components(p);
    p += 3;
  }
  lik_rec[2] -= (foot_min_z-rcap);
  //print_array(lik_rec,config_dim);exit(1);
}

// orients surface (specified by surf_normal)
// is only used externaly, for testing
void ghostmodel::set_surf_rot(const extvec& eas){
  const extvec transl (0,0,0);
  const extvec orient[] = {transl, eas};
  affine_from_orientation(surf_rot, orient);
  //surf_rot.print();
  extvec n (0,0,1);
  surf_rot.mult(n,surf_normal);
  //surf_rot.invert_rigidbody();
  horizontal_flag = false;
}

// rotates torso (by surf_rot that takes surface to horizontal ground)
void ghostmodel::rotate_surf(){
  affine A, B;
  torso_opart->get_frame_A_ground_from_body(A);
  A.get_translation(torso_com);
  set_torso_pos();
  surf_rot.mult(A,B);
  euler_angles_from_affine(B,lik_rec+3);
}

// (x,y,1).dot(plane_zf) produces torso com plane's z component at (x,y) 
void ghostmodel::set_plane_zf(){
  double nx, ny, nz;
  surf_normal.get_components(nx,ny,nz);
  plane_zf.set(-nx/nz,-ny/nz,surf_normal.dot(torso_com)/nz);
}

// transforms points according to rotation 
// of their (vertical) projection onto surface
void ghostmodel::surf_rotate_pos(extvec& pos){
  if(horizontal_flag){return;}
  extvec pos1, pos2, xy1;
  xy1.copy(pos);
  xy1.set_v(2,1);
  pos1.copy(pos);
  pos1.set_v(2,xy1.dot(plane_zf));
  surf_rot.mult(pos1,pos2);
  pos2.subtract(pos1);
  pos.add(pos2);
  //cout<<"pos1: ";pos1.print();cout<<"pos2: ";pos2.print();surf_rot.print();
}

// transforms torso_pos according to torso_com rotation
// note: torso_com is global coord, torso_pos is relative (from free6)
void ghostmodel::surf_rotate_torso_pos(){
  extvec del;
  del.copy(torso_com);
  surf_rotate_pos(del);
  del.subtract(torso_com);
  torso_pos.add(del);
}

// sets torso_pos from torso_com
void ghostmodel::set_torso_pos(){
  torso_pos.copy(torso_com);
  torso_pos.subtract(parent_pos);
  //torso_pos.print();torso_com.print();
  //exit(1);
}

//float th=0;
void ghostmodel::set_surf_rot_from_normal(){
  //if(th<1.){th+=.01;} extvec eas (0,-.2*sin(th),0); set_surf_rot(eas);

  dMatrix3 rot;
  rot_ztov(rot, surf_normal);
  transpose_odematrix(rot);
  surf_rot.set_rotation(rot);

  // to keep torso_com in place 
  extvec tc1; // tc1 = torso_com1
  surf_rot.mult(torso_com,tc1);
  tc1.subtract(torso_com);
  tc1.times(-1);
  surf_rot.translate(tc1);

  set_plane_zf();//plane_zf.print();
  horizontal_flag = false;
}

void ghostmodel::set_gmodel_limb_bends(){
  extvec angles;
  double *p = config+6;
  vector<liklimb*> *limbs, *glimbs;
  limbs = lik->get_limbs();
  glimbs = glik->get_limbs();
  vector<liklimb*>::iterator it, it1;
  it = limbs->begin();
  it1 = glimbs->begin();
  for(;it!=limbs->end();it++){
    angles.set(p);
    p += 3;
    (*it1++)->set_bend((*it)->bend_from_angles(angles));
    //bool bend = (*it)->bend_from_angles(angles);cout << bend << " ";
  }
  //cout<<endl;
}
