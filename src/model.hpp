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
   * The list of neighbours is used for computing the derivatives and is pre-
   * computed during initialization. The neighbours_patch variable has the same
   * role but for a region of the size of the cell patches (all patches have the
   * same size). These variables are computed in Initialize() and do not change
   * at runtime.
   * */
  std::vector<stencil> neighbors, neighbors_patch;

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
  /** Shape parameter: order */
  std::vector<double> S_order;
  /** Shape parameter: angle */
  std::vector<double> S_angle;
  /** Structure tensor */
  std::vector<double> S00, S01;
  /** Q-tensor */
  std::vector<double> Q00, Q01, deltaQ00, deltaQ01, Q00_old, Q01_old;
  /** Nematic field */
  std::vector<double> sumQ00, sumQ01, sumQ00_cnt, sumQ01_cnt;
  /** Elasticity */
  std::vector<double> gam;
  /** Energy penalty for area */
  std::vector<double> mu;
  /** Elongation parameter */
  std::vector<double> delta;

  /** @} */

  /** Domain managment
   * @{ */

  /** Min of the boundaries of the patches and center of mass */
  std::vector<coord> patch_min;
  /** Max of the boundaries of the patches and center of mass */
  std::vector<coord> patch_max;
  /** Size of the patch (set by margin) */
  coord patch_size, patch_margin;
  /** Total number of nodes in a patch */
  unsigned patch_N;
  /** Memory offset for each patch */
  std::vector<coord> offset;
  /** Center-of-mass */
  std::vector<vec<double, 2>> com, com_prev;
  /** Counter to compute com in Fourier space */
  std::vector<std::complex<double>> com_x, com_y;
  /** Precomputed tables for sin and cos (as complex numbers) used in the
   * computation of the com.
   * */
  std::vector<std::complex<double>> com_x_table, com_y_table;

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
  /** shall we perform runtime checks ? */
  bool runtime_check = false;
  /** shall we print some runtime stats ? */
  bool runtime_stats = false;
  /** padding for onscreen output */
  unsigned pad;
  /** name of the inpute file */
  std::string inputname = "";
  /** Switches */
  bool force_delete;
  /** The random number seed */
  unsigned long seed; // picked using a fair die
  /** Flag if seed was set in the arguments, see options.hpp */
  bool set_seed;
  /** Number of predictor-corrector steps */
  unsigned npc = 1;
  /** Relaxation time at initialization */
  unsigned relax_time;
  /** Value of nsubstep to use for initialization */
  unsigned relax_nsubsteps;
  /** Total time spent writing output */
  std::chrono::duration<double> write_duration;
  /** @} */

  /** Simulation parameters
   * @{ */

  /** size of the system */
  vec<unsigned, 2> Size;
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
  /** Margin for the definition of patches */
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
  /** Boudaries for cell generation
   *
   * These are the boundaries (min and max x and y components) of the domain in
   * which the cells are created when the initial config 'random' is choosen.
   * */
  std::vector<unsigned> birth_bdries;
  /** @} */

  /** Cell properties
   * @{ */

  /** Interface thickness */
  double lambda;
  /**  Interaction stength */
  double kappa;
  /** Activity */
  double zeta = 0.;
  /** Internal pressure */
  double P = 0;
  /** Cell-cell friction parameter */
  double f;
  /** Cell-wall friction parameter */
  double f_walls;
  /** Substrate friction parameter */
  double xi;
  /** Adhesion parameter */
  double omega;
  /** Prefered radii (area = pi*R*R) and radius growth */
  std::vector<double> R, dR;
  /** Migration speed */
  double alpha;
  /** Coupling between area and contractility */
  double beta;
  /** Nematic parameters */
  double C, K, D, L;
  /** Pre-computed coefficients */
  double C1, C3;

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
       & auto_name(lambda)
       & auto_name(kappa)
       & auto_name(alpha)
       & auto_name(R)
       & auto_name(xi)
       & auto_name(omega)
       & auto_name(init_config)
       & auto_name(zeta)
       & auto_name(D)
       & auto_name(C)
       & auto_name(K)
       & auto_name(L)
       & auto_name(f)
       & auto_name(f_walls)
       & auto_name(wall_thickness)
       & auto_name(wall_kappa)
       & auto_name(wall_omega)
       & auto_name(walls)
       & auto_name(patch_margin)
       & auto_name(patch_size);
  }

  /** Serialization of parameters (in and out)*/
  template<class Archive>
  void SerializeFrame(Archive& ar)
  {
      ar & auto_name(nphases)
         & auto_name(phi)
         & auto_name(offset)
         & auto_name(area)
         & auto_name(com)
         & auto_name(S_order)
         & auto_name(S_angle)
         & auto_name(Q00)
         & auto_name(Q01)
         & auto_name(velp)
         & auto_name(velf)
         & auto_name(velc)
         & auto_name(patch_min)
         & auto_name(patch_max);
  }

  // =========================================================================
  // Program managment. Implemented in main.cpp

  /** The main loop */
  void Algorithm();

  /** Setup computation */
  void Setup(int, char**);

  /** Do the computation */
  void Run();

  /** Clean after you */
  void Cleanup();

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

  /** Swap two cells in the internal arrays */
  void SwapCells(unsigned n, unsigned m);

  // ===========================================================================
  // Random numbers generation. Implemented in random.cpp

  /** Pseudo random generator */
  std::mt19937 gen;
  //ranlux24 gen;

  /** Return random real, uniform distribution */
  double random_real(double min=0., double max=1.);

  /** Return random real, gaussian distributed */
  double random_normal(double sigma=1.);

  /** Return geometric dist numbers, prob is p */
  unsigned random_geometric(double p);

  /** Return poisson distributed unsigned integers */
  unsigned random_poisson(double lambda);

  /** Return exp distributed unsigned integers */
  int random_exponential(double lambda);

  /** Return random unsigned uniformly distributed */
  unsigned random_unsigned();

  /** Initialize random numbers
   *
   * If CUDA is enabled, alos intialize CUDA random numbers
   * */
  void InitializeRandomNumbers();

  // ==========================================================================
  // Support for cuda. Implemented in cuda.cu

  #ifdef _CUDA_ENABLED

  /** Device(s) propeties
   * @{ */

  /** Obtain (and print) device(s) properties
   *
   * Also checks that the device properties are compatible with the values given
   * in src/cuda.h.
   * */
  void QueryDeviceProperties();

  /** @} */

  /** Pointer to device global memory
   *
   * These pointers reflects the program data strcuture and represents the cor-
   * responding data on the device global memory. All names should be identical
   * to their host counterparts apart from the d_ prefix.
   *
   * @{ */

  double *d_phi, *d_phi_old, *d_V, *d_potential, *d_potential_old, *d_sum,
         *d_square, *d_Q00, *d_Q01, *d_walls,
         *d_walls_laplace, *d_walls_dx, *d_walls_dy, *d_sum_cnt, *d_square_cnt,
         *d_Q00_cnt, *d_Q01_cnt, *d_area, *d_area_cnt, *d_c, *d_S00, *d_S01, *d_S_order,
         *d_S_angle, *d_theta, *d_gam, *d_mu;
  vec<double, 2>  *d_pol, *d_vel, *d_velp, *d_velc, *d_velf, *d_com, *d_com_prev;
  stencil         *d_neighbors, *d_neighbors_patch;
  coord           *d_patch_min, *d_patch_max, *d_offset;
  cuDoubleComplex *d_com_x, *d_com_y, *d_com_x_table, *d_com_y_table;

  /** @} */

  /** Random number generation
   * @{ */

  /** Random states on the device */
  curandState *d_rand_states;

  /** Initialization function */
  void InitializeCuda();

  /** Initialization function for random numbers */
  void InitializeCUDARandomNumbers();

  /** @} */
  /** CUDA device memory managment
    * @{ */

  /** In which direction do we copy data? */
  enum class CopyMemory {
    HostToDevice,
    DeviceToHost
  };

  /** Allocate or free memory? */
  enum class ManageMemory {
    Allocate,
    Free
  };

  /** Implementation for AllocDeviceMemory() and FreeDeviceMemory() */
  void _manage_device_memory(ManageMemory);
  /** Implementation for PutToDevice() and GetFromDevice() */
  void _copy_device_memory(CopyMemory);

  /** Copy data to the device global memory
   *
   * This function is called at the begining of the program just before the main
   * loop but after the system has been initialized.
   * */
  void PutToDevice();

  /** Copy data from the device global memory
   *
   * This function is called every time the results need to be dumped on the disk.
   * */
  void GetFromDevice();

  /** Allocate memory for all device arrays */
  void AllocDeviceMemory();

  /** Allocate memory for all device arrays */
  void FreeDeviceMemory();

  /** @} */

  /** Runtime properties
   * @{ */

  int n_total, n_blocks, n_threads;

  /** @} */

  #endif

  // ===========================================================================
  // Configuration. Implemented in configure.cpp

  /** Add cell with number n at a certain position */
  void AddCell(unsigned n, const coord& center);

  /** Subfunction for AddCell() */
  void AddCellAtNode(unsigned n, unsigned k, const coord& center);

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
  void UpdateAtNode(unsigned, unsigned, bool);

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
  void UpdatePolarization(unsigned, bool);

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

  /** Update the moving patch following each cell */
  void UpdatePatch(unsigned);

  /** Update fields
   *
   * The boolean argument is used to differentiate between the predictor step
   * (true) and subsequent corrector steps.
   * */
  void Update(bool, unsigned=0);

  /** Helper function
   *
   * Update the fields in a square patch that is entirely contained inside the
   * patch, i.e. that is not wrapping around the borders.
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
   * This function is used to updated the fields only in a restricted patch
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

  /** Division flag */
  bool division = false;
  /** Time scale of the division */
  double division_time = 100;
  /** Division rate */
  double division_rate = 0.;
  /** Growth facto before division */
  double division_growth = 1.5;
  /** Relaxation time for division */
  int division_relax_time = 100;
  /** Relaxation time after division */
  int division_refract_time = 100;
  /** Count the number of time steps before the next division */
  std::vector<int> division_counter;

  /** Reset division counter for single cell */
  void ResetDivisionCounter(unsigned);

  /** Make a cell divide
   *
   * TBD
   * */
  void Divide(unsigned i);

  // ===========================================================================
  // Tools

  /** Gives domain coordinate corresponding to a domain index */
  coord GetPosition(unsigned k) const
  { return { GetXPosition(k), GetYPosition(k) }; }

  /** Gives domain x-coordinate corresponding to a domain index */
  unsigned GetXPosition(unsigned k) const
  { return k/Size[1]; }

  /** Gives domain y-coordinate corresponding to a domain index */
  unsigned GetYPosition(unsigned k) const
  { return k%Size[1]; }

  /** Get domain index from domain coordinates */
  unsigned GetIndex(const coord& p) const
  { return p[1] + Size[1]*p[0]; }

  /** Get patch index from domain coordinates */
  unsigned GetPatchIndex(unsigned n, coord p) const
  {
    // get difference to the patch min
    p = (p + Size - patch_min[n])%Size;
    // correct for offset
    p = (p + patch_size - offset[n])%patch_size;
    // remap linearly
    return p[1] + patch_size[1]*p[0];
  }

  /** Get patch index from domain index */
  unsigned GetPatchIndex(unsigned n, unsigned k) const
  {
    return GetPatchIndex(n, GetPosition(k));
  }

  /** Get domain index from patch index */
  unsigned GetIndexFromPatch(unsigned n, unsigned q) const
  {
    // position on the patch
    const coord qpos = { q/patch_size[1], q%patch_size[1] };
    // position on the domain
    const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
    // return domain index
    return GetIndex(dpos);
  }
};

#endif//MODEL_HPP_
