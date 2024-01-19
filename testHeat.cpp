#include "heatProblem.hpp"
#include "inverseforme.hpp"
#include "parameters.hpp"
#include "stringTools.hpp"
#include <iostream>


using namespace std;
using namespace Eigen;



int main(int argc,char* argv[]) {

  std::string paramfile = "params.in";
  Parameters p(paramfile);

  // Check saved configuration
   if (argc > 3)
  {
    p.set("shape", argv[1]);
    p.set("testid", argv[2]);
    p.set("noise", argv[3]);
  }
  else 
  {
    std::cout << "USAGE: ./exec SHAPE TESTID NOISE" << std::endl;
    std::cout << " SHAPE:  1, 2 or 3" << std::endl;
    std::cout << " TESTID: 1, 2 or 3" << std::endl;
    std::cout << " NOISE: a float like 0.02" << std::endl;
    return 1;
  }

  std::string output_dir = p.getString("output dir");
  std::string shape = p.getString("shape");
  std::string testid = p.getString("testid");
  std::string noise = p.getString("noise");
  replace(noise, ".", "_"); // replace dot by underscore

  output_dir += "_shape_" + shape + "_test_" + testid + "_noise_" + noise;
  p.set("output dir", output_dir);

  p.display();
  HeatProblem pb(p);
  
  Inverserr inv(p,pb);
  inv.run_inv_cond_in_any_geom();


 
  
  return 0;
}