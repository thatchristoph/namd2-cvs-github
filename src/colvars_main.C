#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_standalone.h"

int main (int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config_file [restart_prefix]\n";
    return 2;
  }

  colvarproxy_standalone *proxy = 
    new colvarproxy_standalone (std::string (argv[1]),
                                (argc > 2) ? std::string (argv[2]) : std::string (""),
                                std::string ("colvars_out"));

  proxy->colvars->analyse();

  if (proxy->colvars->cv_traj_read_name.size()) {
    // a trajectory was provided
    while (proxy->colvars->read_traj (proxy->colvars->cv_traj_read_name.c_str())) {
      if (proxy->colvars->b_analysis)
        proxy->colvars->analyse();
    }
    if (proxy->colvars->b_analysis)
      proxy->colvars->analyse();
  }
 
  delete proxy;
}
