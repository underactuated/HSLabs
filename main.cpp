#include "core.h"
#include "matrix.h"
#include "visualization.h"
#include "model.h"
#include "player.h"

#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "pergen.h"
#include "periodic.h"

#include "effdata.h" // temp

//using namespace std;
//using namespace rapidxml;

void traverse(xml_node<> *node);

int main(int argc, char *argv[]){

  /*balltreetester tester (10);
  tester.test();
  exit(1);*/

  //srand(time(NULL));
  //srand(4);
  //srand(3);// 4

  /*modelplayer player1;
  player1.load_draw_model("spider.xml");
  //player1.load_draw_model("walker2d5l.xml");
  //player1.load_draw_model("ant.xml");
  exit(1);*/

  modelplayer player0;
  player0.set_play_dt(.02);
  //pergensetup* pgs0 = player0.make_pergensu("pgs_config.txt",1);
  //player0.test(4);exit(1);
  //pergensetup* pgs = player0.make_pergensu("pgs_config.txt",2);
  pergensetup* pgs = player0.make_pergensu("pgs_config.txt",26);
  //pgs->print(2);
  //player0.print_limb_pos0s();exit(1);
  //extvec rec_eas (0,0,-1.571); pgs->set_rec_rotation(rec_eas);
  //player0.test_dynamics(pgs);exit(1);
  player0.speedup_draw(3);//(3);//5);
  player0.set_flag("manual_viewpoint",false);
  //player0.set_flag("texture",true);
  player0.set_flag("smooth_view",true);
  //player0.shift_view(-2,.5,0);
  //player0.shift_view(-1,.5,0);
  //player0.shift_view(0,1,0);
  //player0.set_flag("torso_kicks",true);
  //player0.set_flag("dynamics_from_simulation",true);
  player0.uneven_ground_test();
  //player0.test_lik_solvers();exit(1);
  //player0.play_pergensu(pgs);exit(1);
  //player0.record_per_traj(pgs);exit(1);
  //player0.record_per_traj_sweep(pgs,"step_length",-.5e-6,.5e-6,2);exit(1);
  //player0.record_per_traj_sweep(pgs,"step_length",-.5,.5,1);exit(1);
  //player0.record_per_traj_sweep(pgs,"step_length",-.6,.6,1);exit(1);
  //player0.record_per_traj_sweep(pgs,"step_length",-.5,-.3,1);exit(1);
  //player0.record_per_traj_sweep(pgs,"step_length",-.2,.6,10);exit(1);
  //player0.record_per_traj_sweep(pgs,"period",2,5,7);exit(1);
  //player0.record_pos_control_traj(pgs,0,3.5+10);exit(1);
  //player0.cpc_test(pgs,0);exit(1);
  player0.position_control_test(pgs,0);exit(1);
  //player0.open_loop_test(pgs,0);exit(1);
  //player0.simulate_pergensu(pgs,0);exit(1);
  //string str;player0.pergensu_config_string(pgs,str);cout<<str<<endl;exit(1);
  player0.set_flag("contact_force",true);
  //player0.measure_cot_sweep(pgs,20,"period",13,17,2);exit(1);
  player0.measure_cot_sweep(pgs,20,"period",3,18,15);exit(1);
  //player0.measure_cot_sweep(pgs,20,"step_duration",0,1,10);exit(1);
  //double cot = player0.measure_cot(pgs,20);cout<<"COT = "<<cot<<endl;exit(1);
  periodic per (player0.get_model());
  per.record_trajectory(pgs,20);
  //per.print_trajectory(99,2); exit(1);
  per.compute_dynrecs();
  per.compute_dynrec_ders();
  //per.print_dynrecs(0,0,5);
  player0.get_model()->print();
  //per.print();exit(1);
  //per.print_dynrec(2);exit(1);
  //per.print_dynrecs(3,0,5); //exit(1);
  //per.forcetorque_matrix_sp(2);//exit(1);
  //per.solve_forcetorques(2);exit(1);
  //per.print_dynrecs(3,2,1);
  //per.print_dynrecs(3,12,1);
  //per.print_dynrecs(3);
  per.switch_torso_penalty(1,0+1);
  double work = per.work_over_period(); cout<<"work = "<<work<<endl;exit(1);
  per.check_solve_ft();exit(1);
  //per.print_dynrecs(4);
  //per.print_trajectory(); exit(1);
  //player0.set_play_dt(.02);
  player0.play_pergensu(pgs); 
  delete pgs; exit(1);

  /*kinematicmodel model;
  model.load_fromxml("myant.xml");
  model.recompute_modelnodes();
  periodic per (model);
  per.print();exit(1);*/

  modelplayer player;
  player.load_model("myant.xml");
  //player.load_model("hexapod.xml");
  player.set_flag("manual_viewpoint",false);
  player.test(4);

  /*kinematicmodel model (true);
  //model.load_fromxml("ant.xml");
  model.load_fromxml("myant.xml");
  model.print(); //exit(1);
  //model.recompute_modelnodes();
  //model.test_joints(true);
  model.orient_bodys();
  model.draw();*/

  /*kinematicmodel model (true);
  model.load_fromxml("hexapod.xml");
  model.print(); //exit(1);
  model.orient_bodys();
  model.draw(); exit(1);*/

  return 0;
}

void traverse(xml_node<> *node){
  cout << node->name() << endl;
  xml_node<> *n = node->first_node();
  while(n){
    //cout << n->name() << endl;
    traverse(n);
    n = n->next_sibling();
  }
}


