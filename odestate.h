/////////////////////////////////////////////////
// odestate.h: declaration of structures needed
// for storing, printing and restoring of ODE
// body states, and class for online update of
// configuration time derivatives.
/////////////////////////////////////////////////
#ifndef ODESTATE_H
#define ODESTATE_H

#include "core.h"
#include "model.h"

// Structure odebodydata stores mutable ODE body data.
struct odebodydata{
  dBodyID odebody;
  dVector3 pos, vel, ang_vel; // position, velocity, angular velocity
  dQuaternion quat; // rotation
  odebodydata(dBodyID odebody_){odebody = odebody_;}
  void read(); // from odebody
  void write(); // to odebody
  void print();
};

// Structure odestate stores data of all ODE bodys.
struct odestate{
  list<odebodydata*> obodys;
  odestate(const kinematicmodel* model);
  ~odestate();
  void save(); // from model
  void load(); // to model
  void print();
};

// Class configtimedertrack ( = configuration time derivative
// tracker) keeps track of time derivatives of current
// configuration up to order max_der. Derivatives are
// updated as soon as a new configuration is pushed in.
class configtimedertrack{
  int max_der, config_dim;
  double dt;
  double** ders; // derivatives
  double *new_der, *old_der; // pointers to the last two derivatives
  arrayops ao;
public:
  configtimedertrack(int max_der, int config_dim, double dt);
  configtimedertrack(const configtimedertrack* tdertrack);
  ~configtimedertrack();
  inline double** get_ders(){return ders;}
  void push_config(double* config);
  void print();
  void del_second_der(configtimedertrack* tdertrack, double* del);
private:
  void construct(int max_der, int config_dim, double dt);
  void copy(const configtimedertrack* tdertrack);
  void compute_dern(int n);
};

#endif
