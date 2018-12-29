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

#ifdef DEBUG
#include <fenv.h>
#endif

using namespace std;

/** pure eyecandy */
string title = R"(
                   ______     __          __
                  / ____/__  / /___ _____/ /________
                 / /   / _ \/ / __ `/ __  / ___/ __ \
                / /___/  __/ / /_/ / /_/ / /  / /_/ /
                \____/\___/_/\__,_/\__,_/_/   \____/
                -----------------------------------
                 Celadro: Cells as active droplets
                    (c) 2016-19, Romain Mueller

)";

// =============================================================================

void Model::Algorithm()
{
  // number of steps between two writes
  const unsigned streak_length = nsubsteps*ninfo;

  for(unsigned t=0; t<nsteps; t+=ninfo)
  {
    // write current frame
    if(!no_write and t>=nstart)
    {
      const auto start = chrono::steady_clock::now();

      try
      {
        WriteFrame(t);
      }
      catch(...) {
        cerr << "error" << endl;
        throw;
      }

      write_duration += chrono::steady_clock::now() - start;
    }

    // some verbose
    if(verbose>1) cout << '\n';
    if(verbose)
      cout << "timesteps t = " << setw(pad) << setfill(' ') << right
                                 << t << " to "
                                 << setw(pad) << setfill(' ') << right
                                 << t+ninfo << endl;
    if(verbose>1) cout << string(width, '-') << endl;

    // do the computation
    for(unsigned s=0; s<streak_length; ++s)
    {
      // first sweeps produces estimate of values
      // subsequent sweeps produce corrected values
      for(unsigned i=0; i<=npc; ++i)
        Update(i==0);
    }

#ifdef _CUDA_ENABLED
    // get data from device to host memory
    GetFromDevice();
#endif

    // runtime stats and checks
    try
    {
      if(runtime_stats and verbose>1) RuntimeStats();
      if(runtime_check) RuntimeChecks();
    }
    catch(const error_msg& e)
    {
      throw;
    }
    catch(const warning_msg& e)
    {
      if(stop_at_warning) throw;
      else if(verbose and !no_warning) cerr << "warning: " << e.what() << "\n";
    }
  }

  // finally write final frame
  if(!no_write and nsteps>=nstart) WriteFrame(nsteps);
}

void Model::Setup(int argc, char **argv)
{
  // ========================================
  // Setup

  if(argc<2) throw error_msg("no argument provided. Type -h for help.");
  // parse program options
  ParseProgramOptions(argc, argv);
  // print simulation parameters
  if(verbose)
  {
    cout << "Run parameters" << endl;
    cout << string(width, '=') << endl;
    PrintProgramOptions();
  }

  // ========================================
  // Initialization
  if(verbose) cout << endl << "Initialization" << endl << string(width, '=')
                   << endl;

  // warning and flags
  // no output
  if(no_write and verbose) cout << "warning: output is not enabled." << endl;

  // ----------------------------------------
  // model init
  if(verbose) cout << "model initialization ..." << flush;
  try {
    InitializeRandomNumbers();
    Initialize();
    InitializeNeighbors();
  } catch(...) {
    if(verbose) cout << " error" << endl;
    throw;
  }
  if(verbose) cout << " done" << endl;

  // ----------------------------------------
  // parameters init
  if(verbose) cout << "system initialisation ..." << flush;
  try {
    Configure();
    ConfigureWalls(BC);
  } catch(...) {
    if(verbose) cout << " error" << endl;
    throw;
  }
  if(verbose) cout << " done" << endl;

  // ----------------------------------------
  // multi-threading
  #ifdef _OPENMP
    if(nthreads)
    {
      if(verbose) cout << "multi-threading ... " << flush;
      SetThreads();
      if(verbose) cout << nthreads << " active threads" << endl;
    }
  #endif

  // ----------------------------------------
  // cuda, see random numbers section as well
  #ifdef _CUDA_ENABLED

    if(verbose) cout << "setting up CUDA devices ..." << endl;
    QueryDeviceProperties();
    InitializeCuda();

    if(verbose) cout << "... allocate device memory ...";
    AllocDeviceMemory();
    if(verbose) cout << " done" << endl;

    if(verbose) cout << "... random numbers initialization ..." << flush;
    InitializeCUDARandomNumbers();
    if(verbose) cout << " done" << endl;

    if(verbose) cout << "... copy data to device ...";
    PutToDevice();
    if(verbose) cout << " done" << endl;
  #endif

  // ----------------------------------------
  // write params to file
  if(!no_write)
  {
    if(verbose) cout << "create output directory " << " ...";
    try {
      ClearOutput();
    } catch(...) {
      if(verbose) cout << " error" << endl;
      throw;
    }
    if(verbose) cout << " done" << endl;

    if(verbose)
      cout << "write parameters ...";

    try {
      WriteParams();
    } catch(...) {
      if(verbose) cout << " error" << endl;
      throw;
    }
    if(verbose) cout << " done" << endl;
  }
}

void Model::Run()
{
  // preparation
  if(verbose)   cout << "preparation ... " << flush;
  Pre();
  if(verbose) cout << " done" << endl;

  if(verbose) cout << endl << "Run" << endl << string(width, '=') << "\n\n";

  // print some stats
  PreRunStats();

  // record starting time
  const auto start = chrono::steady_clock::now();
  // run the thing
  Algorithm();
  // record end time
  const auto duration = chrono::steady_clock::now() - start;

  if(verbose) cout << "post-processing ... " << flush;
  Post();
  if(verbose) cout << "done" << endl;

  if(verbose)
  {
    cout << endl << "Statistics" << endl << string(width, '=') << endl;
    cout << "Total run time :                    "
         << chrono::duration_cast<chrono::milliseconds>(duration).count()
            /1000. << " s" << endl;
    cout << "Total time spent writing output :   "
         << chrono::duration_cast<chrono::milliseconds>(write_duration).count()
            /1000. << " s" << endl;
  }
}

void Model::Cleanup()
{
  #ifdef _CUDA_ENABLED
    FreeDeviceMemory();
  #endif
}

/** Program entry */
int main(int argc, char **argv)
{
  // if in debug mode, catch all arithmetic exceptions
  #ifdef DEBUG
    feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif

  // print that beautiful title
  cout << title << endl;

  // do the job
  try {

    Model model;
    model.Setup(argc, argv);
    model.Run();
    model.Cleanup();

  }
  // custom small messages
  catch(const error_msg& e) {
    cerr << argv[0] << ": error: " << e.what() << endl;
    return 1;
  }
  // bad alloc (possibly from memory())
  catch(const bad_alloc& ba) {
    cerr << argv[0] << ": error initializing memory: " << ba.what() << endl;
    return 1;
  }
  // all the rest (mainly from boost)
  catch(const exception& e) {
    cerr << argv[0] << ": " << e.what() << endl;
    return 1;
  }

  return 0;
}
