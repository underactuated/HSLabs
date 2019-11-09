#ifndef ODESTATE_H
#define ODESTATE_H

#include "core.h"
#include "model.h"

struct odebodydata{
  dBodyID odebody;
  dVector3 pos, vel, ang_vel;
  dQuaternion quat;
  odebodydata(dBodyID odebody_){odebody = odebody_;}
  void read(); // from odebody
  void write(); // to odebody
  void print();
};

struct odestate{
  list<odebodydata*> obodys;
  odestate(const kinematicmodel* model);
  ~odestate();
  void save(); // from model
  void load(); // to model
  void print();
};

class configtimedertrack{
  int max_der, config_dim;
  double dt;
  double** ders;
  double *new_der, *old_der;
  arrayops ao;
public:
  configtimedertrack(int max_der, int config_dim, double dt);
  configtimedertrack(const configtimedertrack* tdertrack);
  ~configtimedertrack();
  inline double** get_ders(){return ders;}
  void push_config(double* config);
  void compute_dern(int n);
  void print();
  void del_second_der(configtimedertrack* tdertrack, double* del);
private:
  void construct(int max_der, int config_dim, double dt);
  void copy(const configtimedertrack* tdertrack);
};

#endif
