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

#include "cuda.hpp"

/** Typed malloc on device memory */
template<class T>
inline void tmalloc(T* ptr, size_t length)
{
  cudaMalloc((void**)&ptr, length*sizeof(T));
}

/** Typed memcpy to device memory */
template<class T>
inline void tmemput(T* dest, T*src, size_t length)
{
  cudaMalloc((void**)&ptr, length*sizeof(T));
}

template<class T>
void Model::AllocDeviceMemory()
{
  //
  // must be called after Initialize()
  //

  cmalloc(d_walls, N);
  cmalloc(d_walls_dx, N);
  cmalloc(d_walls_dy, N);
  cmalloc(d_walls_laplace, N);
  cmalloc(d_sum, N);
  cmalloc(d_sum_cnt, N);
  cmalloc(d_square, N);
  cmalloc(d_square_cnt, N);
  cmalloc(d_Px, N);
  cmalloc(d_Py, N);
  cmalloc(d_Theta, N);
  cmalloc(d_Q00, N);
  cmalloc(d_Q01, N);
  cmalloc(d_Px_cnt, N);
  cmalloc(d_Py_cnt, N);
  cmalloc(d_Theta_cnt, N);
  cmalloc(d_Q00_cnt, N);
  cmalloc(d_Q01_cnt, N);
  cmalloc(d_P, N);
  cmalloc(d_P_cnt, N);

  cmalloc(d_phi, nphases*patch_N);
  cmalloc(d_phi_old, nphases*patch_N);
  cmalloc(d_V, nphases*patch_N);
  cmalloc(d_potential, nphases*patch_N);
  cmalloc(d_potential_old, nphases*patch_N);

  cmalloc(d_area, nphases);
  cmalloc(d_area_cnt, nphases);
  cmalloc(d_patch_min, nphases);
  cmalloc(d_patch_max, nphases);
  cmalloc(d_com, nphases);
  cmalloc(d_com_prev, nphases);
  cmalloc(d_pol, nphases);
  cmalloc(d_velp, nphases);
  cmalloc(d_velf, nphases);
  cmalloc(d_velc, nphases);
  cmalloc(d_com_x, nphases);
  cmalloc(d_com_y, nphases);
  cmalloc(d_c, nphases);
  cmalloc(d_S00, nphases);
  cmalloc(d_S01, nphases);
  cmalloc(d_S_order, nphases);
  cmalloc(d_S_angle, nphases);
  cmalloc(d_theta, nphases);
  cmalloc(d_offset, nphases);
}

void Model::FreeDeviceMemory()
{
  cudaFree(d_walls);
  cudaFree(d_walls_dx);
  cudaFree(d_walls_dy);
  cudaFree(d_walls_laplace);
  cudaFree(d_sum);
  cudaFree(d_sum_cnt);
  cudaFree(d_square);
  cudaFree(d_square_cnt);
  cudaFree(d_Px);
  cudaFree(d_Py);
  cudaFree(d_Theta);
  cudaFree(d_Q00);
  cudaFree(d_Q01);
  cudaFree(d_Px_cnt);
  cudaFree(d_Py_cnt);
  cudaFree(d_Theta_cnt);
  cudaFree(d_Q00_cnt);
  cudaFree(d_Q01_cnt);
  cudaFree(d_P);
  cudaFree(d_P_cnt);

  cudaFree(d_phi);
  cudaFree(d_phi_old);
  cudaFree(d_V);
  cudaFree(d_potential);
  cudaFree(d_potential_old);

  cudaFree(d_area);
  cudaFree(d_area_cnt);
  cudaFree(d_patch_min);
  cudaFree(d_patch_max);
  cudaFree(d_com);
  cudaFree(d_com_prev);
  cudaFree(d_pol);
  cudaFree(d_velp);
  cudaFree(d_velf);
  cudaFree(d_velc);
  cudaFree(d_com_x);
  cudaFree(d_com_y);
  cudaFree(d_c);
  cudaFree(d_S00);
  cudaFree(d_S01);
  cudaFree(d_S_order);
  cudaFree(d_S_angle);
  cudaFree(d_theta);
  cudaFree(d_offset);
}

void Model::CopyToDevice()
{
  cudaMemcpy(d_walls, &walls[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_walls_dx, &walls_dx[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_walls_dy, &walls_dy[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_walls_laplace, &walls_laplace[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_sum, &sum[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_sum_cnt, &sum_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_square, &square[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_square_cnt, &square_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Px, &Px[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Py, &Py[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Theta, &Theta[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Q00, &Q00[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Q01, &Q01[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Px_cnt, &Px_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Py_cnt, &Py_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Theta_cnt, &Theta_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Q00_cnt, &Q00_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_Q01_cnt, &Q01_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_P, &P[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_P_cnt, &P_cnt[0], N, cudaMemcpyToDevice);

  for(int i=0; i<patch_N; ++i)
  {
    cudaMemcpy(d_phi+i*patch_N, &phi[i][0], patch_N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_phi_old+i*patch_N, &phi_old[i][0], patch_N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_V+i*patch_N, &V[i][0], patch_N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_potential+i*patch_N, &potential[i][0], patch_N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_potential_old+i*patch_N, &potential_old[i][0], patch_N, cudaMemcpyHostToDevice);
  }

  cudaMemcpy(d_area, &area[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_area_cnt, &area_cnt[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_patch_min, &patch_min[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_patch_max, &patch_max[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_com, &com[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_com_prev, &com_prev[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_pol, &pol[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_velp, &velp[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_velf, &velf[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_velc, &velc[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_com_x, &com_x[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_com_y, &com_y[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_c, &c[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_S00, &S00[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_S01, &S01[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_S_order, &S_order[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_S_angle, &S_angle[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_theta, &theta[0], N, cudaMemcpyToDevice);
  cudaMemcpy(d_offset, &offset[0], N, cudaMemcpyToDevice);
}

void Model::GetFromDevice()
{
}
