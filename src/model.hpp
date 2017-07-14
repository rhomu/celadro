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

#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <vector>
#include <array>
#include "vec.hpp"
#include "stencil.hpp"
#include "serialization.hpp"

/** Type used to represent values on the grid */
using field = std::vector<double>;
/** Grid coordinate */
using coord = vec<unsigned, 2>;

/** Model class
 *
 * This class contains the whole program and is mainly used to be able to
 * scatter the implementation to different files without having to declare
 * every single variable extern. The architecture has not been really thought
 * through and one should not create two of these objects.
 * */
struct Model
{
  /** Simulation variables
   * @{ */

  /** List of neighbours
   *
   * WHY? */
  std::vector<stencil> neighbors;
  /** Phase fields */
  std::vector<field> phi;
  /** Predicted phi in a PC^n step */
  std::vector<field> phi_old;
  /** V = delta F / delta phi */
  std::vector<field> V;
  /** Potential: -V/2 - adv. term */
  std::vector<field> potential;
  /** Predicted potential in a P(C)^n step */
  std::vector<field> potential_old;
  /** Active velocity angle */
  std::vector<double> theta;
  /** Area associated with a phase field */
  std::vector<double> area;
  /** Counter for computing the area */
  std::vector<double> area_cnt;
  /** Sum of phi at each node */
  std::vector<double> sum, sum_cnt;
  /** Sum of square phi at each node */
  std::vector<double> square, square_cnt;
  /** Phase-field for the walls */
  std::vector<double> walls, walls_dx, walls_dy, walls_laplace;
  /** Cell polarity */
  std::vector<vec<double, 2>> pol;
  /** Passive velocities */
  std::vector<vec<double, 2>> velp;
  /** Contraction velocities */
  std::vector<vec<double, 2>> velc;
  /** Friction velocities */
  std::vector<vec<double, 2>> velf;
  /** Total com velocity */
  std::vector<vec<double, 2>> vel;
  /** Center-of-mass */
  std::vector<vec<double, 2>> com, com_prev;
  /** Overall polarization of the tissue */
  std::vector<double> Px, Py, Px_cnt, Py_cnt;
  /** Contractility */
  std::vector<double> c;
  /** Shape parameter: order */
  std::vector<double> S_order;
  /** Shape parameter: angle */
  std::vector<double> S_angle;
  /** Structure tensor */
  std::vector<double> S00, S01;
  /** Polarity tensor */
  std::vector<double> Q00, Q01;
  std::vector<double> Theta;
  std::vector<double> Theta_cnt;
  /** Counters for polarity tensor */
  std::vector<double> Q00_cnt, Q01_cnt;
  /** Internal pressure */
  std::vector<double> P, P_cnt;
   /** Counter to compute com in Fourier space */
  std::vector<std::complex<double>> com_x;
  /** Counter to compute com in Fourier space */
  std::vector<std::complex<double>> com_y;
  /** Precomputed tables for sin and cos (as complex numbers) used in the
   * computation of the com.
   * */
  std::vector<std::complex<double>> com_x_table, com_y_table;
  /** Division flag */
  bool division = false;
  /** Pre-computed coefficients */
  double C1, C2, C3;
  /* Total size of the domain */
  unsigned Size;
  /** @} */

  /** Domain managment
   * @{ */

  /** Memory offset for each domain */
  std::vector<coord> offset;
  /** Min of the boundaries of the domains and center of mass */
  std::vector<coord> domain_min;
  /** Max of the boundaries of the domains and center of mass */
  std::vector<coord> domain_max;

  /** @} */

  /** Program options
   * @{ */

  /** verbosity level
   *
   * 0: no output
   * 1: some output
   * 2: extended output (default)
   * 3: debug
   * */
  unsigned verbose = 2;
  /** compress output? (we use zip) */
  bool compress, compress_full;
  /** name of the run */
  std::string runname;
  /** Output dir (or tmp dir before moving files to the archive) */
  std::string output_dir;
  /** write any output? */
  bool no_write = false;
  /** skip runtime warnings? */
  bool no_warning = false;
  /** are the runtime warnings fatal? (i.e. they do stop the simulation) */
  bool stop_at_warning = false;
  /** padding for onscreen output */
  unsigned pad;
  /** name of the inpute file */
  std::string inputname = "";
  /** Switches */
  bool force_delete;
  /** The random number seed */
  unsigned seed;
  /** Number of predictor-corrector steps */
  unsigned npc = 1;
  /** Relaxation time at initialization */
  unsigned relax_time;
  /** Value of nsubstep to use for initialization */
  unsigned relax_nsubsteps;
  /** Enable tracking? */
  bool tracking = false;
  /** Total time spent writing output */
  std::chrono::duration<double> write_duration;
  /** @} */

  /** Simulation parameters
   * @{ */

  /** size of the system */
  unsigned LX, LY;
  /** total number of nodes */
  unsigned N;
  /** Total number of time steps */
  unsigned nsteps;
  /** Time interval between data outputs */
  unsigned ninfo;
  /** Time at which to start the output */
  unsigned nstart;
  /** number of subdivisions for a time step */
  unsigned nsubsteps;
  /** effective time step */
  double time_step;
  /** Number of phases */
  unsigned nphases;
  /** Boundary conditions flag */
  unsigned BC;
  /** angle in degrees (input variable only) */
  double angle_deg;
  /** Margin for the definition of domains */
  unsigned margin;
  /** Intial configuration */
  std::string init_config;
  /** Noise level for initial configurations */
  double noise = 0;
  /** Wall thickness */
  double wall_thickness = 1.;
  /** Repuslion by the wall */
  double wall_kappa;
  /** Adhesion on the wall */
  double wall_omega;
  /** Division rate */
  double division_rate = 0.;
  /** Relaxation time for division */
  unsigned division_relax_time = 100;
  /** Boudaries for cell generation
   *
   * These are the boundaries (min and max x and y components) of the domain in
   * which the cells are created when the initial config 'random' is choosen.
   * */
  std::vector<unsigned> birth_bdries;
  /** @} */

  /** Cell properties
   * @{ */

  /** Elasticity */
  std::vector<double> gam;
  /** Energy penalty for area */
  std::vector<double> mu;
  /** Interface thickness */
  double lambda;
  /**  Interaction stength */
  double kappa;
  /** Friction parameter */
  double zeta = 0.;
  /** Cell-cell friction parameter */
  double f;
  /** Cell-wall friction parameter */
  double f_walls;
  /** Substrate friction parameter */
  double xi;
  /** Adhesion parameter */
  double omega;
  /** Prefered radius (area = pi*R*R) */
  double R;
  /** Migration speed */
  double alpha;
  /** Coupling between area and contractility */
  double beta;
  /** Parameters for the polarisation dynamics */
  double D=1., J=1.;
  /** Contractility parameters */
  double c0, tauc;

  /** @} */
  /** Multi-threading parameters
   * @{ */
  #ifdef _OPENMP

  /** number of threads */
  unsigned nthreads;

  #endif
  /** @} */

  // ===========================================================================
  // Options. Implemented in options.cpp

  /** the variables map used to collect options */
  opt::variables_map vm;

  /** Set parameters from input */
  void ParseProgramOptions(int ac, char **av);

  /** Output all program options */
  void PrintProgramOptions();

  // ===========================================================================
  // Serialization

  /** Serialization of parameters (in and out)*/
  template<class Archive>
  void SerializeParameters(Archive& ar)
  {
    ar & auto_name(gam)
       & auto_name(mu)
       & auto_name(nphases)
       & auto_name(lambda)
       & auto_name(kappa)
       & auto_name(alpha)
       & auto_name(R)
       & auto_name(xi)
       & auto_name(omega)
       & auto_name(init_config)
       & auto_name(zeta)
       & auto_name(D)
       & auto_name(J)
       & auto_name(f)
       & auto_name(f_walls)
       & auto_name(wall_thickness)
       & auto_name(wall_kappa)
       & auto_name(wall_omega)
       & auto_name(walls);

    ar & auto_name(tracking)
       & auto_name(margin);
  }

  /** Serialization of parameters (in and out)*/
  template<class Archive>
  void SerializeFrame(Archive& ar)
  {
      ar & auto_name(phi)
         & auto_name(area)
         & auto_name(com)
         & auto_name(S_order)
         & auto_name(S_angle)
         & auto_name(pol)
         & auto_name(velp)
         & auto_name(velf)
         & auto_name(velc)
         & auto_name(vel);
      if(tracking) ar
         & auto_name(domain_min)
         & auto_name(domain_max);
  }

  // =========================================================================
  // Program managment. Implemented in main.cpp

  /** The main loop */
  void Algorithm();

  /** Setup computation */
  void Setup(int, char**);

  /** Do the computation */
  void Run();

  // ==========================================================================
  // Writing to file. Implemented in write.cpp

  /** Write current state of the system */
  void WriteFrame(unsigned);

  /** Write run parameters */
  void WriteParams();

  /** Remove old files */
  void ClearOutput();

  /** Create output directory */
  void CreateOutputDir();

  // ==========================================================================
  // Initialization. Implemented in init.cpp

  /** Initialize memory for field */
  void Initialize();

  /** Initialize neighbors list (stencils) */
  void InitializeNeighbors();

  // ===========================================================================
  // Configuration. Implemented in configure.cpp

  /** Add cell with number n at a certain position */
  void AddCell(unsigned n, const coord& center);

  /** Set initial condition for the fields */
  void Configure();

  /** Set initial configuration for the walls */
  void ConfigureWalls();

  // ===========================================================================
  // Run. Implemented in run.cpp

  /** Time step
   *
   * This is the time-stepping function and performs the main computation.
   * */
  void Step();

  /** Prepare before run */
  void Pre();

  /** Prints some stats before running */
  void PreRunStats();

  /** Prints some stats in between ninfo steps */
  void RuntimeStats();

  /** Performs punctual check at runtime */
  void RuntimeChecks();

  /** Post run function */
  void Post();

  /** Subfunction for update */
  void UpdateAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void UpdateFieldsAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void UpdateStructureTensorAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void UpdateFrictionForceAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void SquareAndSumAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void ReinitSquareAndSumAtNode(unsigned);

  /** Compute center of mass of a given phase field */
  void ComputeCoM(unsigned);

  /** Update polarisation of a given field
   *
   * This function updates the polarisation of the cell which give the direction
   * of the active velocity of the cell. We follow reference 10.1101/095133 and
   * define the dynamics as
   *
   *    d theta / dt = J_r torque + 2 D_r eta
   *
   * where eta is gaussian white noise with zero mean and unit variance, see
   * paper for more details. Note that we use the euler-maruyama method, instead
   * of a predictor-corrector method.
   * */
  void UpdatePolarization(unsigned);

  /** Update friction force
   *
   * This function updates the friction force and needs to be in a separate
   * loop because it must be computed after the passive part of the velocity has
   * been computed fully. See paper for more details.
   * */
  void UpdateFrictionForce();

  /** Compute shape parameters
   *
   * This function effectively computes the second moment of area, which ca n be used to
   * fit the shape of a cell to an ellipse.
   * */
  void ComputeShape(unsigned);

  /** Update fields */
  void Update();

  /** Update the window for tracking */
  void UpdateWindow(unsigned);

  /** Helper function
   *
   * Update the fields in a square domain that is entirely contained inside the
   * domain, i.e. that is not wrapping around the borders.
   * */
  template<typename Ret, typename ...Args>
  void UpdateSubDomain(Ret (Model::*)(unsigned, unsigned, Args...),
                       unsigned, unsigned, unsigned, unsigned,
                       unsigned, Args&&... args);
  /** Parallel version */
  template<typename Ret, typename ...Args>
  void UpdateSubDomainP(Ret (Model::*)(unsigned, unsigned, Args...),
                        unsigned, unsigned, unsigned, unsigned,
                        unsigned, Args&&... args);
  /** Helper function
   *
   * This function is used to updated the fields only in a restricted domain
   * around the cell center. One needs to be careful because of the periodic
   * boundary conditions. The template argument is the function used to update
   * the fields at each node (called ***AtNode()).
   * */
  template<typename R, typename ...Args>
  void UpdateDomain(R (Model::*)(unsigned, unsigned, Args...),
                    unsigned, Args&&... args);
  /** Parallel version */
  template<typename R, typename ...Args>
  void UpdateDomainP(R (Model::*)(unsigned, unsigned, Args...),
                     unsigned, Args&&... args);

  // ===========================================================================
  // Division. Implemented in division.cpp

  /** Make a cell divide
   *
   * The strategy for implementing division is to chose a division axis randomly
   * then divide the given cell in two while fixing all the other cells. We then
   * let the two new cells relax while fixing the other cells such that they are
   * able to create a common interface.
   * */
  void Divide(unsigned i);

  // ===========================================================================
  // Tools

  /** Gives x position corresponding to a grid index */
  unsigned GetXPosition(unsigned k) const
  { return k/LY; }

  /** Gives y position corresponding to a grid index */
  unsigned GetYPosition(unsigned k) const
  { return k%LY; }

  /** Get grod index from position */
  unsigned GetDomainIndex(unsigned x, unsigned y) const
  { return y + LY*x; }

  /** Returns total grid size */
  unsigned GetSize() const
  { return LX*LY; }
};

#endif//MODEL_HPP_
