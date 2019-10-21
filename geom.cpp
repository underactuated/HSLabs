#include "core.h"
#include "geom.h"


trimesh::trimesh(int n_vert, int n_tri_, float rgb_[3]){
  n_tri = n_tri_;
  //vertices = new dReal [3*n_vert];
  vertices = new dReal [4*n_vert];
  indices = new dTriIndex [3*n_tri];
  std::copy(rgb_,rgb_+3,rgb);
}

trimesh::~trimesh(){
  delete [] vertices;
  delete [] indices;
}


trimeshmanager::trimeshmanager(dSpaceID space_){
  space = space_;
  clear_buffs();
}

void trimeshmanager::push_vertex(dReal x, dReal y, dReal z){
  dReal a[] = {x,y,z,0};
  for(int i=0;i<4;i++){vertex_buffer.push_back(a[i]);}
  vert_buff_size++;
}

void trimeshmanager::push_triangle(int vi0, int vi1, int vi2){
  int a[] = {vi0,vi1,vi2};
  for(int i=0;i<3;i++){triangle_buffer.push_back(a[i]);}
  tri_buff_size++;
}

void trimeshmanager::set_color(float r, float g, float b){
  rgb[0] = r;
  rgb[1] = g;
  rgb[2] = b;
}

dGeomID trimeshmanager::new_trimesh(){
  trimesh* tm = new trimesh (vert_buff_size, tri_buff_size, rgb);
  dReal *p = tm->vertices;
  list<dReal>::iterator it = vertex_buffer.begin();
  for(int i=0;i<4*vert_buff_size;i++){*p++ = *it++;}
  dTriIndex *p1 = tm->indices;
  list<int>::iterator it1 = triangle_buffer.begin();
  for(int i=0;i<3*tri_buff_size;i++){*p1++ = *it1++;}
  dTriMeshDataID tmdata = dGeomTriMeshDataCreate();
  dGeomTriMeshDataBuildSimple(tmdata, tm->vertices, vert_buff_size, tm->indices, 3*tri_buff_size);
  dGeomID geom = dCreateTriMesh(space,tmdata,0,0,0);
  geom_trimesh[geom] = tm;
  clear_buffs();
  return geom;
}

// destructor is untested
void trimeshmanager::delete_trimesh(dGeomID geom){
  dTriMeshDataID tmdata = dGeomTriMeshGetData(geom);
  dGeomTriMeshDataDestroy (tmdata);
  delete geom_trimesh[geom]; // do we need to destroy geom?
  geom_trimesh.erase(geom);
}

void trimeshmanager::clear_buffs(){
  vertex_buffer.clear();
  triangle_buffer.clear();
  vert_buff_size = 0;
  tri_buff_size = 0;
  set_color(1,1,1);
}

void trimeshmanager::clear(){
  list<dGeomID> geoms;
  map<dGeomID,trimesh*>::iterator it = geom_trimesh.begin();
  for(;it!=geom_trimesh.end();it++){
    geoms.push_back((*it).first);
  }
  list<dGeomID>::iterator it1 = geoms.begin();
  for(;it1!=geoms.end();it1++){delete_trimesh(*it1);}
}

void trimeshmanager::draw(dGeomID geom){
  if(!geom_trimesh.count(geom)){cout<<"ERROR: geom is not registered"<<endl;exit(1);}
  trimesh* tm = geom_trimesh[geom];
  float* col = tm->get_color();
  dsSetColor(col[0],col[1],col[2]);
  int n_tri = tm->n_tri;
  dVector3* vv = new dVector3 [3];
  dReal pos[] = {0,0,0};
  dMatrix3 rot;
  dRSetIdentity(rot);
  for(int i=0;i<n_tri;i++){
    dGeomTriMeshGetTriangle(geom,i,&vv[0],&vv[1],&vv[2]);
    dsDrawTriangle(pos,rot,vv[0],vv[1],vv[2],1);
    //for(int j=0;j<3;j++){print_array<dReal>(vv[j],3);}cout<<endl;
  }
  dsSetColor(1,1,1);
}

void trimeshmanager::print(){
  cout << "trimesh geoms: " << geom_trimesh.size() << endl;
  cout << "vertex buffer: " << vert_buff_size << endl;
  cout << "triangle buffer: " << tri_buff_size << endl;
}


heightfield::heightfield(int nx_, int ny_, double l_){
  nx = nx_;
  ny = ny_;
  l = l_;
  if((nx<1)||(ny<1)||(l<=0)){cout<<"ERROR: wrong params"<<endl;exit(1);}
  heights = new double [(nx+1)*(ny+1)];
  heightcs = new double [nx*ny];
}

heightfield::~heightfield(){
  delete [] heights;
  delete [] heightcs;
}

void heightfield::random_field(double lb, double ub, bool add_flag){
  double *p = heights;
  for(int i=0;i<(nx+1)*(ny+1);i++){
    //*p++ = (lb + double(rand())/RAND_MAX*(ub-lb));
    double h = (lb + double(rand())/RAND_MAX*(ub-lb));
    if(add_flag){*p++ += h;}
    else{*p++ = h;}
  }
  compute_heightcs();
}

dGeomID heightfield::make_geom(double x0_, double y0_, trimeshmanager* trimeshman){
  x0 = x0_;
  y0 = y0_;
  double x, y, z;
  double *p = heights;
  for(int ix=0;ix<=nx;ix++){
    for(int iy=0;iy<=ny;iy++){
      x = (ix-.5*nx)*l + x0;
      y = (iy-.5*ny)*l + y0;
      z = *p++;
      trimeshman->push_vertex(x,y,z);
      //cout<<x<<" "<<y<<" "<<z<<endl;
    }
  }
  p = heightcs;
  for(int ix=0;ix<nx;ix++){
    for(int iy=0;iy<ny;iy++){
      x = (ix-.5*nx+.5)*l + x0;
      y = (iy-.5*ny+.5)*l + y0;
      z = *p++;
      trimeshman->push_vertex(x,y,z);
      //cout<<x<<" "<<y<<" "<<z<<endl;
    }
  }
  int nxy = (nx+1)*(ny+1);
  int iv[4];
  for(int ix=0;ix<nx;ix++){
    for(int iy=0;iy<ny;iy++){
      iv[0] = ix*(ny+1) + iy;
      iv[1] = iv[0] + 1;
      iv[2] = iv[1] + ny + 1;
      iv[3] = iv[2] - 1;
      int ic = nxy + ix*ny + iy;
      for(int i=0;i<4;i++){
	trimeshman->push_triangle(ic,iv[(i+1)%4],iv[i]);
      }
    }
  }
  trimeshman->set_color(.7,.4,.1);
  //trimeshman->print();
  return trimeshman->new_trimesh();
}

void heightfield::compute_heightcs(){
  double *p = heights, *p1 = heights + nx*(ny+1);
  for(int i=0;i<=ny;i++){*p++ = 0; *p1++ = 0;}
  p = heights;
  for(int i=0;i<=nx;i++){*p = 0; p += ny; *p++ = 0;}
  for(int ix=0;ix<nx;ix++){
    for(int iy=0;iy<ny;iy++){
      double hc = 0;
      for(int i=0;i<4;i++){
	hc += *hpointer(ix+int(i/2),iy+(i%2),0);
      }
      *hpointer(ix,iy,1) = hc/4;
    }
  }
}

double* heightfield::hpointer(int ix, int iy, bool c_flag){
  int n = ny + 1;
  double* p = heights;
  if(c_flag){n--; p = heightcs;}
  return p + ix*n + iy;
}

void heightfield::slope_field(double h0, double h1){
  for(int ix=0;ix<=nx;ix++){
    double h = h0 + (h1-h0)*double(ix)/nx;
    for(int iy=0;iy<=ny;iy++){
      *hpointer(ix,iy,0) = h;
    }
  }
  compute_heightcs();
}

double heightfield::get_h(double x, double y){
  int ix = -1, iy = -1;
  double dx, dy;
  nearby_cell_coords(x,y,ix,iy,dx,dy);
  if(ix<0 || ix>=nx || iy<0 || iy>=ny){return 0;}

  int iq = -1;
  if (dy>=0 && dy>=fabs(dx)) {iq = 0;}
  else if (dx>=0 && dx>=fabs(dy)) {iq = 1;}
  else if (dy<=0 && -dy>=fabs(dx)) {iq = 2;}
  else if (dx<=0 && -dx>=fabs(dy)) {iq = 3;}
  //cout<<x<<" "<<y<<" "<<ix<<" "<<iy<<" "<<iq<<" "<<dx<<" "<<dy<<endl;

  double z01c[3];
  ixyq_to_z01c(ix,iy,iq,z01c);
  double z0 = z01c[0], z1 = z01c[1], zc = z01c[2];
  double s0 = (z0-z1)/l, s1 = (z0+z1-2*zc)/l;
  const int f0[] = {-1,0,1,0};
  const int f1[] = {0,1,0,-1};

  double h = zc + dx*(s0*f0[iq]+s1*f1[iq]) + dy*(s0*f1[iq]-s1*f0[iq]);
  return h;
}

void heightfield::nearby_cell_coords(double x, double y, int& ix, int& iy, double& dx, double& dy){
  x -= x0;
  y -= y0;
  ix = floor(x/l+.5*nx);
  iy = floor(y/l+.5*ny);
  double xc = (ix-.5*nx+.5)*l;
  double yc = (iy-.5*ny+.5)*l;
  dx = x-xc;
  dy = y-yc;
}

void heightfield::ixyq_to_z01c(int ix, int iy, int iq, double* z01c){
  int ix0 = -1, iy0 = -1, ix1 = -1, iy1 = -1;
  switch (iq) {
  case 0: ix0 = ix; iy0 = iy+1; ix1 = ix+1; iy1 = iy+1; break;
  case 1: ix0 = ix+1; iy0 = iy+1; ix1 = ix+1; iy1 = iy; break;
  case 2: ix0 = ix+1; iy0 = iy; ix1 = ix; iy1 = iy; break;
  case 3: ix0 = ix; iy0 = iy; ix1 = ix; iy1 = iy+1; break;
  }
  //cout<<"ix0="<<ix0<<" iy0="<<iy0<<" ix1="<<ix1<<" iy1="<<iy1<<endl;
  z01c[0] = *hpointer(ix0,iy0,false);
  z01c[1] = *hpointer(ix1,iy1,false);
  z01c[2] = *hpointer(ix,iy,true);
}

void heightfield::ripple_field(double h0, double h1, bool add_flag){
  for(int ix=0;ix<=nx;ix++){
    double h = h0 + (h1-h0)*(ix % 2);
    for(int iy=0;iy<=ny;iy++){
      double* p = hpointer(ix,iy,0);
      if (add_flag) {*p += h;} else {*p = h;}
      //*hpointer(ix,iy,0) = h;
    }
  }
  compute_heightcs();
}

void heightfield::ridge_field(double h0, double h1){
  for(int ix=0;ix<=nx;ix++){
    double h = h0 + (h1-h0)*double(nx-abs(2*ix-nx))/nx;
    for(int iy=0;iy<=ny;iy++){
      *hpointer(ix,iy,0) = h;
    }
  }
  compute_heightcs();
}

void heightfield::tan_field(double phi, double f){
  for(int ix=0;ix<=nx;ix++){
    double h = f*tan(phi*double(ix)/nx);
    for(int iy=0;iy<=ny;iy++){
      *hpointer(ix,iy,0) = h;
    }
  }
  compute_heightcs();
}

void heightfield::tanh_field(double h1){
  for(int ix=0;ix<=nx;ix++){
    double h = h1*(tanh(3*double(2*ix-nx)/nx)+1)/2;
    for(int iy=0;iy<=ny;iy++){
      *hpointer(ix,iy,0) = h;
    }
  }
  compute_heightcs();
}

void heightfield::gauss_field(double h1){
  for(int ix=0;ix<=nx;ix++){
    double a = double(2*ix-nx)/nx;
    double h = h1*exp(-2.5*a*a);
    for(int iy=0;iy<=ny;iy++){
      *hpointer(ix,iy,0) = h;
    }
  }
  compute_heightcs();
}
