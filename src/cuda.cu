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
#include "cuda.h"
#include "model.hpp"

using namespace std;

__global__
static void seed_rand(curandState *state, unsigned long seed)
{
  int id = threadIdx.x;
  curand_init(seed, id, 0, &state[id]);
}
 
/** Typed malloc/free on device memory
  *
  * We implement malloc and free in the same function to avoid copy-pasting
  * mistakes.
  */
template<class T>
void malloc_or_free(T*& ptr, size_t len, Model::ManageMemory which)
{
  if(which==Model::ManageMemory::Allocate)
    cudaMalloc((void**)&ptr, len*sizeof(T));
  else
    cudaFree(ptr);
}

/** Typed memcpy to device memory
  *
  * We implement both directions in the same function to avoid copy-pasting
  * mistakes.
  */
template<class T, class U>
void bidirectional_memcpy(T* device,
                          U* host,
                          size_t len,
                          Model::CopyMemory dir)
{
  if(dir==Model::CopyMemory::HostToDevice)
    cudaMemcpy(device, (U*)host, len*sizeof(T), cudaMemcpyHostToDevice);
  else
    cudaMemcpy((U*)host, device, len*sizeof(T), cudaMemcpyDeviceToHost);
}

void Model::_manage_device_memory(ManageMemory which)
{
  malloc_or_free(d_walls, N, which);
  malloc_or_free(d_walls_dx, N, which);
  malloc_or_free(d_walls_dy, N, which);
  malloc_or_free(d_walls_laplace, N, which);
  malloc_or_free(d_sum, N, which);
  malloc_or_free(d_sum_cnt, N, which);
  malloc_or_free(d_square, N, which);
  malloc_or_free(d_square_cnt, N, which);
  malloc_or_free(d_Px, N, which);
  malloc_or_free(d_Py, N, which);
  malloc_or_free(d_Theta, N, which);
  malloc_or_free(d_Q00, N, which);
  malloc_or_free(d_Q01, N, which);
  malloc_or_free(d_Px_cnt, N, which);
  malloc_or_free(d_Py_cnt, N, which);
  malloc_or_free(d_Theta_cnt, N, which);
  malloc_or_free(d_Q00_cnt, N, which);
  malloc_or_free(d_Q01_cnt, N, which);
  malloc_or_free(d_P, N, which);
  malloc_or_free(d_P_cnt, N, which);
  malloc_or_free(d_phi, nphases*patch_N, which);
  malloc_or_free(d_phi_old, nphases*patch_N, which);
  malloc_or_free(d_V, nphases*patch_N, which);
  malloc_or_free(d_potential, nphases*patch_N, which);
  malloc_or_free(d_potential_old, nphases*patch_N, which);

  malloc_or_free(d_area, nphases, which);
  malloc_or_free(d_area_cnt, nphases, which);
  malloc_or_free(d_patch_min, nphases, which);
  malloc_or_free(d_patch_max, nphases, which);
  malloc_or_free(d_com, nphases, which);
  malloc_or_free(d_com_prev, nphases, which);
  malloc_or_free(d_pol, nphases, which);
  malloc_or_free(d_velp, nphases, which);
  malloc_or_free(d_velf, nphases, which);
  malloc_or_free(d_velc, nphases, which);
  malloc_or_free(d_com_x, nphases, which);
  malloc_or_free(d_com_y, nphases, which);
  malloc_or_free(d_c, nphases, which);
  malloc_or_free(d_S00, nphases, which);
  malloc_or_free(d_S01, nphases, which);
  malloc_or_free(d_S_order, nphases, which);
  malloc_or_free(d_S_angle, nphases, which);
  malloc_or_free(d_theta, nphases, which);
  malloc_or_free(d_offset, nphases, which);
  malloc_or_free(d_gam, nphases, which);
  malloc_or_free(d_mu, nphases, which);

  malloc_or_free(d_com_x_table, Size[0], which);
  malloc_or_free(d_com_y_table, Size[0], which);

  // random number generation
  malloc_or_free(d_rand_states, N, which);
}

void Model::_copy_device_memory(CopyMemory dir)
{
  bidirectional_memcpy(d_walls, &walls[0], N, dir);
  bidirectional_memcpy(d_walls_dx, &walls_dx[0], N, dir);
  bidirectional_memcpy(d_walls_dy, &walls_dy[0], N, dir);
  bidirectional_memcpy(d_walls_laplace, &walls_laplace[0], N, dir);
  bidirectional_memcpy(d_sum, &sum[0], N, dir);
  bidirectional_memcpy(d_sum_cnt, &sum_cnt[0], N, dir);
  bidirectional_memcpy(d_square, &square[0], N, dir);
  bidirectional_memcpy(d_square_cnt, &square_cnt[0], N, dir);
  bidirectional_memcpy(d_Px, &Px[0], N, dir);
  bidirectional_memcpy(d_Py, &Py[0], N, dir);
  bidirectional_memcpy(d_Theta, &Theta[0], N, dir);
  bidirectional_memcpy(d_Q00, &Q00[0], N, dir);
  bidirectional_memcpy(d_Q01, &Q01[0], N, dir);
  bidirectional_memcpy(d_Px_cnt, &Px_cnt[0], N, dir);
  bidirectional_memcpy(d_Py_cnt, &Py_cnt[0], N, dir);
  bidirectional_memcpy(d_Theta_cnt, &Theta_cnt[0], N, dir);
  bidirectional_memcpy(d_Q00_cnt, &Q00_cnt[0], N, dir);
  bidirectional_memcpy(d_Q01_cnt, &Q01_cnt[0], N, dir);
  bidirectional_memcpy(d_P, &P[0], N, dir);
  bidirectional_memcpy(d_P_cnt, &P_cnt[0], N, dir);

  for(unsigned i=0; i<nphases; ++i)
  {
    bidirectional_memcpy(d_phi+i*patch_N, &phi[i][0], patch_N, dir);
    bidirectional_memcpy(d_phi_old+i*patch_N, &phi_old[i][0], patch_N, dir);
    bidirectional_memcpy(d_V+i*patch_N, &V[i][0], patch_N, dir);
    bidirectional_memcpy(d_potential+i*patch_N, &potential[i][0], patch_N, dir);
    bidirectional_memcpy(d_potential_old+i*patch_N, &potential_old[i][0], patch_N, dir);
  }

  bidirectional_memcpy(d_area, &area[0], nphases, dir);
  bidirectional_memcpy(d_area_cnt, &area_cnt[0], nphases, dir);
  bidirectional_memcpy(d_patch_min, &patch_min[0], nphases, dir);
  bidirectional_memcpy(d_patch_max, &patch_max[0], nphases, dir);
  bidirectional_memcpy(d_com, &com[0], nphases, dir);
  bidirectional_memcpy(d_com_prev, &com_prev[0], nphases, dir);
  bidirectional_memcpy(d_pol, &pol[0], nphases, dir);
  bidirectional_memcpy(d_velp, &velp[0], nphases, dir);
  bidirectional_memcpy(d_velf, &velf[0], nphases, dir);
  bidirectional_memcpy(d_velc, &velc[0], nphases, dir);
  bidirectional_memcpy(d_com_x, &com_x[0], nphases, dir);
  bidirectional_memcpy(d_com_y, &com_y[0], nphases, dir);
  bidirectional_memcpy(d_c, &c[0], nphases, dir);
  bidirectional_memcpy(d_S00, &S00[0], nphases, dir);
  bidirectional_memcpy(d_S01, &S01[0], nphases, dir);
  bidirectional_memcpy(d_S_order, &S_order[0], nphases, dir);
  bidirectional_memcpy(d_S_angle, &S_angle[0], nphases, dir);
  bidirectional_memcpy(d_theta, &theta[0], nphases, dir);
  bidirectional_memcpy(d_offset, &offset[0], nphases, dir);
  bidirectional_memcpy(d_gam, &gam[0], nphases, dir);
  bidirectional_memcpy(d_mu, &mu[0], nphases, dir);

  bidirectional_memcpy(d_com_x_table, &com_x_table[0], Size[0], dir);
  bidirectional_memcpy(d_com_y_table, &com_y_table[0], Size[1], dir);
}

void Model::AllocDeviceMemory()
{
  _manage_device_memory(ManageMemory::Allocate);
}

void Model::FreeDeviceMemory()
{
  _manage_device_memory(ManageMemory::Free);
}

void Model::PutToDevice()
{
  _copy_device_memory(CopyMemory::HostToDevice);
}

void Model::GetFromDevice()
{
  _copy_device_memory(CopyMemory::DeviceToHost);
}

void Model::QueryDeviceProperties()
{
  int devCount;
  cudaGetDeviceCount(&devCount);

  if(devCount>1) throw error_msg("multiple Cuda devices not supported.");
  if(devCount==0) throw error_msg("no cuda device found.");

  cudaDeviceProp DeviceProperties;
  cudaGetDeviceProperties(&DeviceProperties, 0);

  if(verbose)
  {
    const int kb = 1024;
    const int mb = kb * kb;

    cout << "... device " << DeviceProperties.name 
         << " (" << DeviceProperties.major << "." << DeviceProperties.minor << ")" << endl;
    cout << "    ... global memory:        " << DeviceProperties.totalGlobalMem / mb 
         << "mb" << endl;
    cout << "    ... shared memory:        " << DeviceProperties.sharedMemPerBlock / kb
         << "kb" << endl;
    cout << "    ... constant memory:      " << DeviceProperties.totalConstMem / kb
         << "kb" << endl;
    cout << "    ... block registers:      " << DeviceProperties.regsPerBlock << endl;
    cout << "    ... warp size:            " << DeviceProperties.warpSize << endl;
    cout << "    ... threads per block:    " << DeviceProperties.maxThreadsPerBlock << endl;
    cout << "    ... max block dimensions: [ " << DeviceProperties.maxThreadsDim[0]
         << ", " << DeviceProperties.maxThreadsDim[1]  << ", " 
         << DeviceProperties.maxThreadsDim[2] << " ]" << endl;
    cout << "    ... max grid dimensions:  [ " << DeviceProperties.maxGridSize[0]
         << ", " << DeviceProperties.maxGridSize[1]  << ", " << DeviceProperties.maxGridSize[2] 
         << " ]" << endl;
  }

  // some checks
  if(DeviceProperties.warpSize>WarpSize)
    throw error_msg("warp size is incompatble with device value. "
                    "See src/cuda.h.");

  if(DeviceProperties.warpSize<WarpSize)
    throw warning_msg("warp size does not match with device value. "
                      "See src/cuda.h.");

  if(DeviceProperties.maxThreadsPerBlock>ThreadsPerBlock)
    throw error_msg("number of threads per block is incompatible with device value. "
                    "See src/cuda.h.");

  if(DeviceProperties.maxThreadsPerBlock<ThreadsPerBlock)
    throw error_msg("number of threads per block does not match with device value. "
                    "See src/cuda.h.");
}

void Model::InitializeCUDARandomNumbers()
{
  seed_rand<<<1, n_total>>>(d_rand_states, seed);
}

void Model::InitializeCuda()
{
  n_total   = static_cast<int>(nphases*patch_N);
  n_blocks  = 1 + (n_total-1)/ThreadsPerBlock;
  n_threads = ThreadsPerBlock;
}
