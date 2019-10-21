using namespace Eigen;

typedef SparseMatrix<double> SpMat; // column-major sparse matrix

struct dynpart{
  int id, parent_id;
  extvec com_pos, joint_pos, foot_pos;
  odepart* opart;
  modelnode* mnode;
  double mass;
  affine A_inertia; // check if affine should be replaced with 3x3 matrices
  bool foot_flag;
  dynpart(odepart* opart_){opart = opart_;}
  void setup(map<modelnode*,int>& mnode_id_map);
  void print();
  extvec* get_com_pos(){return &com_pos;}
  extvec* get_joint_pos(){return &joint_pos;}
  extvec* get_foot_pos(){return &foot_pos;}
  void recompute();
  affine* get_A_ground();
  affine* get_inertia_tensor();
  int get_parent_id(){return parent_id;}
  double get_mass(){return mass;}
  void setup_foot(set<modelnode*>& foot_set);
  bool if_foot(){return foot_flag;}
  void get_joint_zaxis(extvec& axis);
private:
  void set_joint_pos();
  void set_com_pos();
  void set_inertial_params();
  void set_foot_pos();
};

// maybe replace it with more compact/efficient represention
struct dynrecord{
  int n; // number of parts
  int nf; // mumber of feet
  extvec *pos, *jpos, *vel, *mom, *mom_rate, *acc, *ust, *ang_vel, *ang_mom, *ang_mom_rate, *fpos, *jzaxis; // ust stands for u*sin(theta) where u is the unit axis vector
  affine *rot;
  bool* contacts;
  dynrecord(int n, int nf);
  ~dynrecord();
  void initialize(vector<dynpart*>& dynparts, double rcap);
  void print(int id);
  void compute_ders(int stage, dynrecord* prev_rec, dynrecord* next_rec, double dt, periodic* per);
  extvec* get_mom_rate(){return mom_rate;}
  extvec* get_ang_mom_rate(){return ang_mom_rate;}
  extvec* get_pos(){return pos;}
  extvec* get_jpos(){return jpos;}
  extvec* get_jzaxis(){return jzaxis;}
  void set_forcetorque_system(SpMat& B, VectorXd& f, int* parentis, double* masses);
  int get_ncontacts();
  void set_forcetorque_system_contacts(SpMat& B, int* footis);
  void set_forcetorque_system_contacts(SpMat& B, int* footis, bool contact_feet_flag);
private:
  void compute_ders(extvec *p_der, extvec *p_func_prev, extvec *p_func_next, double dt);
  void compute_ang_mom(periodic* per);
  void compute_mom(double* masses);
  void ftsys_forces(int i, int pi, SpMat& B, VectorXd& f);
  void ftsys_torques(int i, int pi, SpMat& B, VectorXd& f);
  void ftsys_gravity(int i, VectorXd& f, double* ms);
  void ftsys_contact_forces(int i, int ci, SpMat& B);
  void ftsys_contact_torques(int fi, int i, int ci, SpMat& B);
};
