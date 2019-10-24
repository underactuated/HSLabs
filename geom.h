#ifndef GEOM_H
#define GEOM_H

#include <list>
#include <map>
#include <ode/ode.h>

#include "core.h"

struct trimesh{
  dReal* vertices;
  dTriIndex* indices;
  int n_tri;
  float rgb[3];
  trimesh(int n_vert, int n_tri, float rgb[3]);
  trimesh(){}
  ~trimesh();
  inline float* get_color(){return rgb;}
};

class trimeshmanager{
  dSpaceID space;
  list<dReal> vertex_buffer;
  list<int> triangle_buffer;
  int vert_buff_size, tri_buff_size;
  float rgb[3];
  map<dGeomID,trimesh*> geom_trimesh;
public:
  trimeshmanager(dSpaceID space);
  void push_vertex(dReal x, dReal y, dReal z);
  void push_triangle(int vi0, int vi1, int vi2);
  void set_color(float r, float g, float b);
  void disable();
  dGeomID new_trimesh();
  void delete_trimesh(dGeomID geom);
  void draw(dGeomID geom);
  void clear();
  void print();
private:
  void clear_buffs();
};

class heightfield{
  int nx, ny;	// size of field in l units
  double l;	// length of (square) cell side
  double *heights, *heightcs;	// y-wise ordering
  double x0, y0;
public:
  heightfield(int nx, int ny, double l);
  ~heightfield();
  void random_field(double lb, double ub, bool add_flag);
  dGeomID make_geom(double x, double y, trimeshmanager* trimeshman);
  void slope_field(double h0, double h1);
  double get_h(double x, double y);
  void ripple_field(double h0, double h1, bool add_flag);
  void ridge_field(double h0, double h1);
  void tan_field(double phi, double f);
  void tanh_field(double h1);
  void gauss_field(double h1);
private:
  void compute_heightcs();
  double* hpointer(int ix, int iy, bool c_flag);
  void nearby_cell_coords(double x, double y, int& ix, int& iy, double& dx, double& dy);
  void ixyq_to_z01c(int ix, int iy, int iq, double* z01c);
};

#endif
