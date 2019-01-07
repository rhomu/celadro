/*
 * This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "header.hpp"
#include "model.hpp"

using namespace std;

void Model::ParseProgramOptions(int ac, char **av)
{
  // ===========================================================================
  // Declaration

  // options allowed only in the command line
  opt::options_description generic("Generic options");
  generic.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v", opt::value<unsigned>(&verbose)->implicit_value(2),
     "verbosity level (0=none, 1=little, 2=normal, 3=debug)")
    ("input,i", opt::value<string>(&inputname),
     "input directory")
    ("force-delete,f", opt::bool_switch(&force_delete),
     "force deletion of existing output files")
#ifdef _OPENMP
    ("threads,t",
     opt::value<unsigned>(&nthreads)->default_value(0)->implicit_value(1),
     "number of threads (0=no multithreading, 1=OpenMP default, "
     ">1=your favorite number)")
#endif
    ("compress,c", opt::bool_switch(&compress),
     "compress individual files using zip")
    ("no-write", opt::bool_switch(&no_write),
     "disable file output (for testing purposes)")
    ("no-warning", opt::bool_switch(&no_warning),
     "disable model specific runtime warnings")
    ("stop-at-warning", opt::bool_switch(&stop_at_warning),
     "runtime warnings interrupt the algorithm")
    ("check", opt::bool_switch(&runtime_check),
     "perform runtime checks")
    ("stat", opt::bool_switch(&runtime_stats),
     "print runtime stats");

  // model specific options
  opt::options_description simulation("Simulation options");
  simulation.add_options()
    ("LX", opt::value<unsigned>(&Size[0]),
     "# of nodes in the x direction")
    ("LY", opt::value<unsigned>(&Size[1]),
     "# of nodes in the y direction")
    ("seed", opt::value<unsigned long>(&seed),
     "set seed for random number generation (random if unset)")
    ("nstart", opt::value<unsigned>(&nstart)->default_value(0u),
     "time at which to start the output")
    ("bc", opt::value<unsigned>(&BC)->default_value(0u),
     "boundary conditions flag (0=pbc, 1=box, 2=channel, 3=ellipse)")
    ("nsteps", opt::value<unsigned>(&nsteps),
     "iterate this many steps in total")
    ("nsubsteps", opt::value<unsigned>(&nsubsteps)->default_value(1u),
     "subdivision of a single time step "
     "(effective time step is 1/nsubsteps)")
    ("ninfo", opt::value<unsigned>(&ninfo),
     "save frame every so many steps")
    ("lambda", opt::value<double>(&lambda),
      "Interface thickness parameter")
    ("kappa", opt::value<double>(&kappa),
      "Interaction strength")
    ("npc", opt::value<unsigned>(&npc)->default_value(1u),
      "Number of predictor-corrector steps")
    ("margin", opt::value<unsigned>(&margin)->default_value(0u),
      "Margin for the definition of restricted domains (if 0: update full box)")
    ("omega", opt::value<double>(&omega),
      "Adhesion parameter")
    ("wall-thickness", opt::value<double>(&wall_thickness),
      "Wall thickness (typical decay length)")
    ("wall-kappa", opt::value<double>(&wall_kappa)->default_value(kappa),
      "Wall repulsion")
    ("wall-omega", opt::value<double>(&wall_omega)->default_value(0.),
      "Wall adhesion")
    ("align-polarization-to", opt::value<int>(&align_polarization_to),
     "Align polarization to velocity (=0) or pressure force (=1)")
    ("align-nematic-to", opt::value<int>(&align_nematic_to),
     "Align nematic tensor to velocity (=0), pressure force (=1), or shape (=2)");

  // cells options
  opt::options_description cells("Cells options");
  cells.add_options()
    ("nphases", opt::value<unsigned>(),
      "Number of phases")
    ("gamma", opt::value<double>(),
     "Elastic constant of each phase (array)")
    ("mu", opt::value<double>(),
      "Energy penalty for area of each phase (array)")
    ("R", opt::value<double>(),
      "Preferred radius (defines area Pi*R*R)")
    ("xi", opt::value<double>(),
      "Substrate friction parameter")
    ("alpha", opt::value<double>()->default_value(0),
     "Strength of propulsion")
    ("J-pol", opt::value<double>()->default_value(0),
     "Nematic flow alignment strength")
    ("S-pol", opt::value<double>()->default_value(0),
      "Norm of the polarisation vector")
    ("D-pol", opt::value<double>()->default_value(0),
      "Polarisation noise strength")
    ("K-pol", opt::value<double>()->default_value(0),
     "elastic constant for the polarisation")
    ("zetaQ", opt::value<double>()->default_value(0),
     "Activity from internal nematic tensor")
    ("zetaS", opt::value<double>()->default_value(0),
     "Activity from shape")
    ("K-nem", opt::value<double>()->default_value(0),
     "elastic constant for the nematic")
    ("J-nem", opt::value<double>()->default_value(0),
     "Nematic flow alignment strength")
    ("W-nem", opt::value<double>()->default_value(0),
     "Strength of vorticity torque")
    ("D-nem", opt::value<double>()->default_value(0),
      "Nematic noise strength")
    ("S-nem", opt::value<double>()->default_value(0),
      "Order of the nematic tensors")
    ("division", opt::value<bool>()->default_value(false),
      "Does the cell and its descendents divide?")
    ("division-rate", opt::value<double>()->default_value(0),
      "Rate of division (in time steps)")
    ("division-time", opt::value<double>()->default_value(100),
      "Time scale for the division");

  // init config options
  opt::options_description init("Initial configuration options");
  init.add_options()
    ("config", opt::value<string>(&init_config),
      "Initial configuration")
    ("relax-time", opt::value<unsigned>(&relax_time)->default_value(0u),
      "Relaxation time steps at initialization.")
    ("noise", opt::value<double>(&noise),
      "Noise level for initial nematic angle, in (0,1).")
    ("cross-ratio", opt::value<double>(&cross_ratio),
      "Ratio of the size of the cross compared to the domain size (for BC=4)")
    ("wound-ratio", opt::value<double>(&wound_ratio),
      "Ratio of the size of the wound open space compared to the domain size (for BC=5)")
    ("tumor-ratio", opt::value<double>(&tumor_ratio),
      "Ratio of the size of the tumor compared to the domain size (for BC=6)")
    ("birth-boundaries", opt::value<vector<unsigned>>(&birth_bdries)->multitoken(),
     "Boundaries in which the cells are created "
     "when the initial configuration 'random' is choosed. "
     "Format: {min x, max, x, min y, max y}")
    ("relax-nsubsteps", opt::value<unsigned>(&relax_nsubsteps)->default_value(0u),
      "Value of nsubsteps to use at initial relaxation (0 means use nsubsteps).");

  // ===========================================================================
  // Parsing

  // command line options
  opt::options_description cmdline_options;
  cmdline_options.add(generic);
  // config file options
  opt::options_description config_file_options;
  config_file_options.add(simulation).add(init);
  // cells file options
  opt::options_description cells_file_options;
  cells_file_options.add(cells);

  // first unnamed argument is the input directory
  opt::positional_options_description p;
  p.add("input", 1);

  // reintialize vm in case we run this function twice
  vm = opt::variables_map();

  // parse first the cmd line
  opt::store(
    opt::command_line_parser(ac, av)
    .options(cmdline_options)
    .positional(p)
    .run(), vm);
  opt::notify(vm);

  // print help msg and exit
  if(vm.count("help"))
  {
    cout << cmdline_options << endl;
    cout << config_file_options << endl;
    cout << cells_file_options << endl;
    exit(0);
  }

  // parse main config file
  if(inputname.empty() or not fs::is_directory(inputname))
    throw error_msg("please provide an input directory. Type -h for help.");
  else
  {
    string config_file_name = inputname + "/config.dat";
    std::fstream file(config_file_name.c_str(), std::fstream::in);
    if(!file.good())
      throw error_msg("can not open runcard file ", config_file_name);
    opt::store(opt::parse_config_file(file, config_file_options), vm);
    opt::notify(vm);
  }

  // parse cells files
  for(auto& p: fs::directory_iterator(inputname + "/cells/"))
  {
    opt::variables_map vm;

    std::fstream file(p.path().c_str(), std::fstream::in);
    if(!file.good())
      throw error_msg("can not open runcard file ", p.path());
    opt::store(opt::parse_config_file(file, cells_file_options), vm);

    // add cells properties
    nphases += vm["nphases"].as<unsigned>();

    gam.resize(nphases, vm["gamma"].as<double>());
    mu.resize(nphases, vm["mu"].as<double>());
    R.resize(nphases, vm["R"].as<double>());
    target_R.resize(nphases, vm["R"].as<double>());
    xi.resize(nphases, vm["xi"].as<double>());
    alpha.resize(nphases, vm["alpha"].as<double>());
    Dpol.resize(nphases, vm["D-pol"].as<double>());
    Spol.resize(nphases, vm["S-pol"].as<double>());
    Jpol.resize(nphases, vm["J-pol"].as<double>());
    Kpol.resize(nphases, vm["K-pol"].as<double>());
    Dnem.resize(nphases, vm["D-nem"].as<double>());
    Snem.resize(nphases, vm["S-nem"].as<double>());
    Jnem.resize(nphases, vm["J-nem"].as<double>());
    Knem.resize(nphases, vm["K-nem"].as<double>());
    Wnem.resize(nphases, vm["W-nem"].as<double>());
    zetaQ.resize(nphases, vm["zetaQ"].as<double>());
    zetaS.resize(nphases, vm["zetaS"].as<double>());
    division.resize(nphases, vm["division"].as<bool>());
    division_time.resize(nphases, vm["division-time"].as<double>());
    division_rate.resize(nphases, vm["division-rate"].as<double>());
    types.resize(nphases, p.path().stem());

    // store variable map in case we need it later
    vm_cells.emplace_back(p.path().stem(), vm);
  }

  if(nphases==0)
    throw error_msg("no cells found, add the relevant config files in the "
                    "cell/ directory");

  // ===========================================================================
  // Fixing some secondary values

  // init random numbers?
  set_seed = vm.count("seed");

  // compute the correct padding
  pad = inline_str(nsteps).length();

  // compute effective time step
  time_step = 1./nsubsteps;

  // set nstart to the next correct frame (round above)
  if(nstart%ninfo) nstart = (1u+nstart/ninfo)*ninfo;
}

/** Print variables from variables_map
  *
  * from: https://gist.github.com/gesquive/8673796
  */
void print_vm(const opt::variables_map& vm, unsigned padding)
{
  for (opt::variables_map::const_iterator it = vm.begin(); it != vm.end(); ++it)
  {
    // pass if defaulted or empty
    if (vm[it->first].defaulted() || it->second.defaulted()) continue;
    if (((boost::any)it->second.value()).empty()) continue;

    std::cout << std::left << std::setw(floor(padding/2)) << it->first;

    /*if (((boost::any)it->second.value()).empty()) {
      std::cout << "(empty)";
    }
    if (vm[it->first].defaulted() || it->second.defaulted()) {
      std::cout << "(default)";
    }*/

    std::cout << std::right << std::setw(ceil(padding/2));

    bool is_char;
    try {
      boost::any_cast<const char*>(it->second.value());
      is_char = true;
    } catch (const boost::bad_any_cast &) {
      is_char = false;
    }
    bool is_str;
    try {
      boost::any_cast<std::string>(it->second.value());
      is_str = true;
    } catch (const boost::bad_any_cast &) {
      is_str = false;
    }

    if (((boost::any)it->second.value()).type() == typeid(int)) {
      std::cout << vm[it->first].as<int>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(unsigned)) {
      std::cout << vm[it->first].as<unsigned>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(size_t)) {
      std::cout << vm[it->first].as<size_t>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
      std::cout << (vm[it->first].as<bool>() ? "true" : "false") << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(double)) {
      std::cout << vm[it->first].as<double>() << std::endl;
    } else if (((boost::any)it->second.value()).type()
               == typeid(vector<double>)) {
      std::cout << vec2str(vm[it->first].as<vector<double>>()) << std::endl;
    } else if (((boost::any)it->second.value()).type()
               == typeid(vector<unsigned>)) {
      std::cout << vec2str(vm[it->first].as<vector<unsigned>>()) << std::endl;
    } else if (is_char) {
      std::cout << vm[it->first].as<const char *>() << std::endl;
    } else if (is_str) {
      std::string temp = vm[it->first].as<std::string>();
      if (temp.size()) {
        std::cout << temp << std::endl;
      } else {
        std::cout << "true" << std::endl;
      }
    } else { // Assumes that the only remainder is vector<string>
      try {
        auto vect = vm[it->first].as<std::vector<std::string> >();
        uint i = 0;
        for (auto oit=vect.begin();
            oit != vect.end(); oit++, ++i) {
          std::cout << "\r> " << it->first
                    << "[" << i << "]=" << (*oit) << std::endl;
        }
      } catch (const boost::bad_any_cast &) {
        std::cout << "UnknownType("
                  << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
      }
    }
  }
}

void Model::PrintProgramOptions()
{
  print_vm(vm, width);

  for(const auto& p: vm_cells)
  {
    cout << "\nCell type " << p.first << ":" << endl;
    print_vm(p.second, width);
  }
}
