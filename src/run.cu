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

/*
 * CUDA implementation: this part of the program is highly experimental ! Use
 * with love and care.
 *
 * Some remarks:
 *   * The implementation reflects the single/multi-core implementation in
 *     src/run.cpp for ease of future development. Hence it is not very
 *     sophisticated and can be probably improved by a substantial margin.
 *     Note however that the two implementations are independent and that any
 *     change made in src/run.cpp must be implemented here independently.
 *   * Some device specific properties can be tuned in src/cuda.h.
 *
 * Implementation details:
 *   * Most of the code can be copy-pasted from the corresponding functions in
 *     run.cpp. In the device memory, all arrays are stored contiguously and we
 *     generally denote by m=0,...,nphases*patch_N this total index. Other
 *     indices are then obtained from it.
 */

#include "header.hpp"
#include "model.hpp"
#include "derivatives.hpp"
#include "cuda.h"
#include "reduce.h"

using namespace std;

// =============================================================================
// kernels

__global__
void UpdateFieldsCuda(
    double *phi,
    double *V,
    vec<double, 2> *pol,
    vec<double, 2> *velp,
    vec<double, 2> *velc,
    vec<double, 2> *velf,
    coord Size,
    coord patch_size,
    coord patch_margin,
    unsigned patch_N,
    coord *patch_min,
    coord *patch_max,
    coord *offset,
    double *sum,
    double *square,
    double *Q00,
    double *Q01,
    double *P,
    double *Px,
    double *Py,
    stencil *neighbors,
    stencil *neighbors_patch,
    double *walls,
    double *walls_laplace,
    double *walls_dx,
    double *walls_dy,
    double alpha,
    double xi,
    double C1,
    double C3,
    double f,
    double f_walls,
    double zeta,
    double kappa,
    double omega,
    double wall_kappa,
    double wall_omega,
    int nphases,
    int n_total
    )
{
  // ---------------------------------------------------------------------------
  // build original indices

  const int m = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned n = static_cast<unsigned>(m)/patch_N;
  const unsigned q = static_cast<unsigned>(m)%patch_N;

  const coord qpos = { q/patch_size[1], q%patch_size[1] };
  const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
  const auto  k    = dpos[1] + Size[1]*dpos[0];

  const auto& s    = neighbors[k];
  const auto& sq   = neighbors_patch[q];

  // ---------------------------------------------------------------------------

  // cell properties
  const auto& p  = phi[m];
  const auto  ll = laplacian(&phi[n*patch_N], sq);
  const auto  dx = derivX(&phi[n*patch_N], sq);
  const auto  dy = derivY(&phi[n*patch_N], sq);
  // all-cells properties
  const auto  ls   = laplacian(sum, s);
  const auto  dxs  = derivX(sum, s);
  const auto  dys  = derivY(sum, s);
  const auto  dxp0 = derivX(Px, s);
  const auto  dyp0 = derivY(Px, s);
  const auto  dxp1 = derivX(Py, s);
  const auto  dyp1 = derivY(Py, s);

  const auto  lw = walls_laplace[k];
  const auto dxw = walls_dx[k];
  const auto dyw = walls_dy[k];

  // delta F / delta phi
  const double force = (
      + C1*(
        // repulsion term
        + kappa*p*(square[k]-p*p)
        + wall_kappa*p*walls[k]*walls[k]
        )
      - C3*(
        // adhesion term
        + omega*(ls-ll)
        + wall_omega*lw
        )
      );

  // potential
  V[m]     = force;

  atomicAdd(&velp[n][0], dx*force);
  atomicAdd(&velp[n][1], dy*force);
  // contractility force
  atomicAdd(&velc[n][0], (P[k]+zeta*Q00[k])*dx + zeta*Q01[k]*dy );
  atomicAdd(&velc[n][1], zeta*Q01[k]*dx + (P[k]-zeta*Q00[k])*dy );
  // friction force
  atomicAdd(&velf[n][0], f*alpha/xi*pol[n][0]*(dx*(dxs-dx)+dy*(dys-dy))
    - f*alpha/xi*(dx*(dxp0-pol[n][0]*dx)+dy*(dyp0-pol[n][0]*dy))
    //+ f_walls*alpha/xi*pol[n][0]*(dx*dxw+dy*dyw);
    - dyw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw));
  atomicAdd(&velf[n][1], f*alpha/xi*pol[n][1]*(dx*(dxs-dx)+dy*(dys-dy))
    - f*alpha/xi*(dx*(dxp1-pol[n][1]*dx)+dy*(dyp1-pol[n][1]*dy))
    //+ f_walls*alpha/xi*pol[n][1]*(dx*dxw+dy*dyw);
    + dxw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw));
}

__global__
void UpdateCuda(
    double *phi,
    double *phi_old,
    double *V,
    double *potential,
    double *potential_old,
    double *gam,
    double *mu,
    double *area,
    double *area_cnt,
    double *theta,
    vec<double, 2> *pol,
    vec<double, 2> *vel,
    vec<double, 2> *velp,
    vec<double, 2> *velc,
    vec<double, 2> *velf,
    coord *patch_min,
    coord *patch_max,
    coord *offset,
    coord Size,
    coord patch_size,
    coord patch_margin,
    vec<double, 2> *com,
    cuDoubleComplex *com_x,
    cuDoubleComplex *com_y,
    double *sum,
    double *square,
    double *P,
    double *Theta,
    double *Q00,
    double *Q01,
    double *Px,
    double *Py,
    double *sum_cnt,
    double *square_cnt,
    double *P_cnt,
    double *Theta_cnt,
    double *Q00_cnt,
    double *Q01_cnt,
    double *Px_cnt,
    double *Py_cnt,
    cuDoubleComplex *com_x_table,
    cuDoubleComplex *com_y_table,
    stencil *neighbors,
    stencil *neighbors_patch,
    double alpha,
    double time_step,
    double xi,
    double C1,
    double C2,
    double J,
    double D,
    unsigned nphases,
    unsigned patch_N,
    bool store,
    int n_total,
    curandState *rand_states
    )
{
  // ---------------------------------------------------------------------------
  // build original indices

  const int m = blockIdx.x*blockDim.x + threadIdx.x;
  if(m>=n_total) return;

  const unsigned n = static_cast<unsigned>(m)/patch_N;
  const unsigned q = static_cast<unsigned>(m)%patch_N;

  // GetIndexFromPatch:
  const coord qpos = { q/patch_size[1], q%patch_size[1] };
  const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
  const auto  k    = dpos[1] + Size[1]*dpos[0];

  const auto& sq   = neighbors_patch[q];

  // ---------------------------------------------------------------------------
  // Update

  {
    const auto p  = phi[m];
    const auto a  = area[n];
    const auto dx = derivX(&phi[n*patch_N], sq);
    const auto dy = derivY(&phi[n*patch_N], sq);
    const auto ll = laplacian(&phi[n*patch_N], sq);

    potential[m] = (
      // free energy term
      -.5*V[m]
      -.5*(
        + C1*gam[n]*p*(1.-p)*(1.-2.*p)
        - 2.*mu[n]*(1.-a/C2)*2.*p
        - 2.*gam[n]*ll
      )
      // advection term
      - (alpha*pol[n][0]+vel[n][0])*dx/xi
      - (alpha*pol[n][1]+vel[n][1])*dy/xi
      );
  }

  // store values
  if(store)
  {
    potential_old[m] = potential[m];
    phi_old[m]       = phi[m];
  }

  // predictor-corrector
  {
    const double p = phi_old[m]
                     + time_step*.5*(potential[m] + potential_old[m]);

    // update for next call
    phi[m]       = p;

    atomicAdd(&sum_cnt[k], p);
    atomicAdd(&area_cnt[n], p*p);
    atomicAdd(&square_cnt[k], p*p);

    const auto cp = make_cuDoubleComplex(p, 0.);
    const auto rx = cuCmul(com_x_table[dpos[0]], cp);
    const auto ry = cuCmul(com_y_table[dpos[1]], cp);
    atomicAdd(&com_x[n].x, rx.x); 
    atomicAdd(&com_x[n].y, rx.y); 
    atomicAdd(&com_y[n].x, ry.x); 
    atomicAdd(&com_y[n].y, ry.y); 
  }

  // save total velocity
  if(q==0)
  {
    vel[n] = velp[n] + velf[n] + velc[n];
  }
}

__global__
void Update2Cuda(
    double *phi,
    double *phi_old,
    double *V,
    double *potential,
    double *potential_old,
    double *gam,
    double *mu,
    double *area,
    double *area_cnt,
    double *theta,
    vec<double, 2> *pol,
    vec<double, 2> *vel,
    vec<double, 2> *velp,
    vec<double, 2> *velc,
    vec<double, 2> *velf,
    coord *patch_min,
    coord *patch_max,
    coord *offset,
    coord Size,
    coord patch_size,
    coord patch_margin,
    vec<double, 2> *com,
    cuDoubleComplex *com_x,
    cuDoubleComplex *com_y,
    double *sum,
    double *square,
    double *P,
    double *Theta,
    double *Q00,
    double *Q01,
    double *Px,
    double *Py,
    double *sum_cnt,
    double *square_cnt,
    double *P_cnt,
    double *Theta_cnt,
    double *Q00_cnt,
    double *Q01_cnt,
    double *Px_cnt,
    double *Py_cnt,
    cuDoubleComplex *com_x_table,
    cuDoubleComplex *com_y_table,
    stencil *neighbors,
    stencil *neighbors_patch,
    double alpha,
    double time_step,
    double xi,
    double C1,
    double C2,
    double J,
    double D,
    unsigned nphases,
    unsigned patch_N,
    bool store,
    int n_total,
    curandState *rand_states
    )
{
  // ---------------------------------------------------------------------------
  // build original indices

  const int m = blockIdx.x*blockDim.x + threadIdx.x;
  if(m>=n_total) return;

  const unsigned n = static_cast<unsigned>(m)/patch_N;
  const unsigned q = static_cast<unsigned>(m)%patch_N;

  // GetIndexFromPatch:
  const coord qpos = { q/patch_size[1], q%patch_size[1] };
  const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
  const auto  k    = dpos[1] + Size[1]*dpos[0];

  const auto& sq   = neighbors_patch[q];

  // -------------------------------------------------------------------------
  // UpdatePolarization

  if(q==0)
  {
    // -------------------------------------------------------------------------
    // UpdatePolarization

    if(store)
    {
      double v[2] = { vel[n][0], vel[n][1] };

      // ...the norm of the passive velocity
      const double ni = sqrt(abs(v[0]*v[0]+v[1]*v[1]));

      // alignement torque
      const double torque = -ni*atan2(v[0]*pol[n][1]-v[1]*pol[n][0],
                             v[0]*pol[n][0]+v[1]*pol[n][1]);

      // ...euler-marijuana update
      theta[n] += time_step*J*torque + sqrt(time_step)*D
                  *curand_normal(&rand_states[n]);

      // update polarisation and contractility
      pol[n] = { cos(theta[n]), sin(theta[n]) };
    }

    // -------------------------------------------------------------------------
    // ComputeCoM

    {
      const auto s = make_cuDoubleComplex(double(Size[0]*Size[1]), 0.);
      const auto tx = cuCdiv(com_x[n], s);
      const auto mx = atan2(tx.y, tx.x) + Pi;
      const auto ty = cuCdiv(com_y[n], s);
      const auto my = atan2(ty.y, ty.x) + Pi;
      com[n] = { mx/2./Pi*Size[0], my/2./Pi*Size[1] };
    }

    // -------------------------------------------------------------------------
    // ComputeShape

    // -------------------------------------------------------------------------
    // UpdatePatch

    // obtain the new location of the patch min and max
    const coord new_min = {
      (static_cast<unsigned>(round(com[n][0])) + Size[0] - patch_margin[0])%Size[0],
      (static_cast<unsigned>(round(com[n][1])) + Size[1] - patch_margin[1])%Size[1]
    };
    const coord new_max = {
      (static_cast<unsigned>(round(com[n][0])) + patch_margin[0] - 1u)%Size[0],
      (static_cast<unsigned>(round(com[n][1])) + patch_margin[1] - 1u)%Size[1]
    };
    coord displacement = (Size + new_min - patch_min[n])%Size;
    // I guess there is somehthing better than this...
    if(displacement[0]==Size[0]-1u) displacement[0] = patch_size[0]-1u;
    if(displacement[1]==Size[1]-1u) displacement[1] = patch_size[1]-1u;

    // update offset and patch location
    offset[n]    = ( offset[n] + patch_size - displacement )%patch_size;
    patch_min[n] = new_min;
    patch_max[n] = new_max;

    // reset values
    com_x[n] = com_y[n] = make_cuDoubleComplex(0., 0.);
    area[n] = 0.;
  }

  // reinit values
  sum[k]    = 0;
  square[k] = 0;
  P[k]      = 0;
  Theta[k]  = 0;
  Q00[k]    = 0;
  Q01[k]    = 0;
  Px[k]     = 0;
  Py[k]     = 0;

  if(q==0)
  {
    velp[n] = {0., 0.};
    velc[n] = {0., 0.};
    velf[n] = {0., 0.};
  }

  // -------------------------------------------------------------------------
  // 

  {
    //const double p = phi[m];
    //atomicAdd(&P_cnt[k]     , p*c[n]);
    /*atomicAdd(&Q00_cnt[k]   , p*.5*(pol[n][0]*pol[n][0]-pol[n][1]*pol[n][1]));
    atomicAdd(&Q01_cnt[k]   , p*pol[n][0]*pol[n][1]);
    atomicAdd(&Px_cnt[k]    , p*pol[n][0]);
    atomicAdd(&Py_cnt[k]    , p*pol[n][1]);
    atomicAdd(&Theta_cnt[k] , p*theta[n]);*/
  }
}

template<class T>
__global__
void TestReduce(T *in, T *out, int N)
{
  T sum = { 0 };

  for(int i = blockIdx.x * blockDim.x + threadIdx.x;
      i < N; i += blockDim.x * gridDim.x)
    sum += in[i];

  blockReduceSum(sum);

  if (threadIdx.x==0) out[blockIdx.x]=sum;
}

// =============================================================================
// Update function

void Model::Update(bool store)
{
/*
  {
    int n_total  = 1024;
    int n_blocks = 1+(n_total-1)/n_threads;
    //int n_blocks = min(65535, 1+(n_total-1)/n_threads);
    int n_stride = 1+(n_total-1)/(n_blocks*n_threads);

    using T = vec<int, 2>;

    T *gam = new T[n_total];
    T *sum = new T[n_blocks];
    T *d_gam;
    T *d_sum;
    cudaMalloc(&d_gam, n_total*sizeof(T));
    cudaMalloc(&d_sum, n_blocks*sizeof(T));

    for(int i=0; i<n_total; ++i) gam[i] = { 1 };//{ 1, 2 , 1, 2};

    cout << "\nTotal: " << n_total << " Blocks: " << n_blocks << " Stride: " << n_stride << "\n";

    cudaMemcpy(d_gam, gam, n_total*sizeof(T), cudaMemcpyHostToDevice);

    TestReduce<<<n_blocks, n_threads>>>(d_gam, d_sum, n_total);
    TestReduce<<<1, n_threads>>>(d_sum, d_sum, n_blocks);

    cudaMemcpy(gam, d_gam, n_total*sizeof(T), cudaMemcpyDeviceToHost);
    cudaMemcpy(sum, d_sum, n_blocks*sizeof(T), cudaMemcpyDeviceToHost);

    //for(int i=0; i<nphases; ++i) cout << gam[i] << ' ';
    cout << "Total: " << sum[0] << '\n';

    cudaFree(d_gam);
    cudaFree(d_sum);
  }

  exit(0);
*/

  UpdateFieldsCuda<<<n_blocks, n_threads>>>(
    d_phi,
    d_V,
    d_pol,
    d_velp,
    d_velc,
    d_velf,
    Size,
    patch_size,
    patch_margin,
    patch_N,
    d_patch_min,
    d_patch_max,
    d_offset,
    d_sum,
    d_square,
    d_Q00,
    d_Q01,
    d_P,
    d_Px,
    d_Py,
    d_neighbors,
    d_neighbors_patch,
    d_walls,
    d_walls_laplace,
    d_walls_dx,
    d_walls_dy,
    alpha,
    xi,
    C1,
    C3,
    f,
    f_walls,
    zeta,
    kappa,
    omega,
    wall_kappa,
    wall_omega,
    nphases,
    n_total);

  UpdateCuda<<<n_blocks, n_threads>>>(d_phi,
    d_phi_old,
    d_V,
    d_potential,
    d_potential_old,
    d_gam,
    d_mu,
    d_area,
    d_area_cnt,
    d_theta,
    d_pol,
    d_vel,
    d_velp,
    d_velc,
    d_velf,
    d_patch_min,
    d_patch_max,
    d_offset,
    Size,
    patch_size,
    patch_margin,
    d_com,
    d_com_x,
    d_com_y,
    d_sum,
    d_square,
    d_P,
    d_Theta,
    d_Q00,
    d_Q01,
    d_Px,
    d_Py,
    d_sum_cnt,
    d_square_cnt,
    d_P_cnt,
    d_Theta_cnt,
    d_Q00_cnt,
    d_Q01_cnt,
    d_Px_cnt,
    d_Py_cnt,
    d_com_x_table,
    d_com_y_table,
    d_neighbors,
    d_neighbors_patch,
    alpha,
    time_step,
    xi,
    C1,
    C2,
    J,
    D,
    nphases,
    patch_N,
    store,
    n_total,
    d_rand_states
    );

  Update2Cuda<<<n_blocks, n_threads>>>(d_phi,
    d_phi_old,
    d_V,
    d_potential,
    d_potential_old,
    d_gam,
    d_mu,
    d_area,
    d_area_cnt,
    d_theta,
    d_pol,
    d_vel,
    d_velp,
    d_velc,
    d_velf,
    d_patch_min,
    d_patch_max,
    d_offset,
    Size,
    patch_size,
    patch_margin,
    d_com,
    d_com_x,
    d_com_y,
    d_sum,
    d_square,
    d_P,
    d_Theta,
    d_Q00,
    d_Q01,
    d_Px,
    d_Py,
    d_sum_cnt,
    d_square_cnt,
    d_P_cnt,
    d_Theta_cnt,
    d_Q00_cnt,
    d_Q01_cnt,
    d_Px_cnt,
    d_Py_cnt,
    d_com_x_table,
    d_com_y_table,
    d_neighbors,
    d_neighbors_patch,
    alpha,
    time_step,
    xi,
    C1,
    C2,
    J,
    D,
    nphases,
    patch_N,
    store,
    n_total,
    d_rand_states
    );

  swap(d_sum, d_sum_cnt);
  swap(d_area, d_area_cnt);
  swap(d_square, d_square_cnt);

/*
  swap(d_P, d_P_cnt);
  swap(d_Theta, d_Theta_cnt);
  swap(d_Q00, d_Q00_cnt);
  swap(d_Q01, d_Q01_cnt);
  swap(d_Px, d_Px_cnt);
  swap(d_Py, d_Py_cnt);*/
}
