struct odebody{
  dBodyID body;
  dVector3 pos, vel, ang_vel;
  dQuaternion quat;
  odebody(dBodyID body_){body = body_;}
  void read(); // from body
  void write(); // to body
  void print();
};

struct odestate{
  list<odebody*> obodys;
  odestate(kinematicmodel* model);
  ~odestate();
  void save(); // from model
  void load(); // to model
  void print();
};

struct configtimedertrack{
  int max_der, config_dim;
  double dt;
  double** ders;
  double *new_der, *old_der;
  arrayops ao;
  configtimedertrack(int max_der, int config_dim, double dt);
  configtimedertrack(configtimedertrack* tdertrack);
  ~configtimedertrack();
  void push_config(double* config);
  void compute_dern(int n);
  void print();
  void del_second_der(configtimedertrack* tdertrack, double* del);
  double** get_ders(){return ders;}
private:
  void construct(int max_der, int config_dim, double dt);
  void copy(configtimedertrack* tdertrack);
};


