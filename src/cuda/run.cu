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
 *   * Some device specific properties can be tuned in src/cuda/cuda.cuh.
 */

#include "cuda/reduce.cuh"

__global__
void UpdateFieldsAtNode(
    double **phi,
    double **V,
    vec<double, 2> *pol,
    vec<double, 2> *velp,
    vec<double, 2> *velc,
    vec<double, 2> *velf,
    coord Size,
    coord patch_size,
    coord patch_margin,
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
    double wall_omega
    int nphases
    )
{
  const int ntot   = blockIdx.x*blockDim.x + threadIdx.x;
  if(ntot>nphases*patch_N) return;

  const unsigned n = static_cast<unsigned>(ntot)/patch_N;
  const unsigned q = static_cast<unsigned>(ntot)%patch_N;

  const coord qpos = { q/patch_size[1], q%patch_size[1] };
  const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
  const auto  k    = dpos[1] + Size[1]*dpos[0];

  const auto& s    = neighbors[k];
  const auto& sq   = neighbors_patch[q];

  // ---------------------------------------------------------------------------

  // cell properties
  const auto& p  = phi[ntot];
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
  V[ntot]     = force;

  vec<double, 6> sums = {
    dx*force,
    dy*force,
    (P[k]+zeta*Q00[k])*dx + zeta*Q01[k]*dy,
    zeta*Q01[k]*dx + (P[k]-zeta*Q00[k])*dy,
    f*alpha/xi*pol[n][0]*(dx*(dxs-dx)+dy*(dys-dy))
      - f*alpha/xi*(dx*(dxp0-pol[n][0]*dx)+dy*(dyp0-pol[n][0]*dy))
      - dyw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw),
    f*alpha/xi*pol[n][1]*(dx*(dxs-dx)+dy*(dys-dy))
      - f*alpha/xi*(dx*(dxp1-pol[n][1]*dx)+dy*(dyp1-pol[n][1]*dy))
      + dxw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw)
  };

  sums = blockReduceSum(sums);

  if(q==0)
  {
    velp[n] = { sums[0], sums[1] };
    velc[n] = { sums[2], sums[3] };
    velf[n] = { sums[4], sums[5] };
  }
}

__global__
void UpdateAtNode(
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
    complex<double> *com_x,
    complex<double> *com_y,
    double *sum,
    double *square,
    double *P,
    double *Theta,
    double *Q00,
    double *Q01,
    double *Px,
    double *Py,
    complex<double> *com_x_table,
    complex<double> *com_y_table,
    stencil *neighbors_patch,
    double alpha,
    double time_step,
    double xi,
    double C1,
    double C2,
    double J,
    double D,
    bool store
    )
{
  const int ntot   = blockIdx.x*blockDim.x + threadIdx.x;
  if(ntot>nphases*patch_N) return;

  const unsigned n = static_cast<unsigned>(ntot)/patch_N;
  const unsigned q = static_cast<unsigned>(ntot)%patch_N;

  const coord qpos = { q/patch_size[1], q%patch_size[1] };
  const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
  const auto  k    = dpos[1] + Size[1]*dpos[0];

  const auto& s    = neighbors[k];
  const auto& sq   = neighbors_patch[q];

  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Update

  {
    const auto  p  = phi[ntot];
    const auto  a  = area[n];
    const auto  dx = derivX(&phi[n*patch_N], sq);
    const auto  dy = derivY(&phi[n*patch_N], sq);
    const auto  ll = laplacian(&phi[n*patch_N], sq);

    potential[n][q] = (
      // free energy term
      -.5*V[ntot]
      -.5*(
        + C1*gam[n]*p*(1.-p)*(1.-2.*p)
        - 2.*mu[n]*(1.-a/C2)*2.*p
        - 2.*gam[n]*ll
      )
      // advection term
      - (alpha*pol[n][0]+velp[n][0]+velc[n][0]+velf[n][0])*dx/xi
      - (alpha*pol[n][1]+velp[n][1]+velc[n][1]+velf[n][1])*dy/xi
      );
  }

  // store values
  if(store)
  {
    potential_old[ntot] = potential[ntot];
    phi_old[ntot]       = phi[ntot];
  }

  // predictor-corrector
  {
    double p = phi_old[ntot]
               + time_step*.5*(potential[ntot] + potential_old[ntot]);

    // update for next call
    phi[ntot]    = p;
    area_cnt[n] += p*p;
    com_x[n]    += com_x_table[dpos[0]]*p;
    com_y[n]    += com_y_table[dpos[1]]*p;
  }

  // -------------------------------------------------------------------------
  // UpdatePolarization

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
    // -------------------------------------------------------------------------
    // UpdatePolarization

    if(store)
    {
      // Polarization...
      array<double, 2> v = {
        velp[n][0] + velf[n][0] + velc[n][0],
        velp[n][1] + velf[n][1] + velc[n][1]
      };

      // ...the norm of the passive velocity
      const double ni = sqrt(abs(v[0]*v[0]+v[1]*v[1]));

      // alignement torque
      double torque = -ni*atan2(v[0]*pol[n][1]-v[1]*pol[n][0],
          v[0]*pol[n][0]+v[1]*pol[n][1]);

      // ...euler-marijuana update
      theta[n] += time_step*J*torque + sqrt(time_step)*D*random_normal();

      // update polarisation and contractility
      pol[n] = { cos(theta[n]), sin(theta[n]) };
    }

    // -------------------------------------------------------------------------
    // ComputeCoM

    const auto mx = arg(com_x[n]/static_cast<double>(Size[0]*Size[1])) + Pi;
    const auto my = arg(com_y[n]/static_cast<double>(Size[0]*Size[1])) + Pi;
    com[n] = { mx/2./Pi*Size[0], my/2./Pi*Size[1] };

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
    com_x[n] = com_y[n] = area[n] = 0.;
  }
}

void Model::Update()
{
  int ntot = static_cast<int>(nphases*patch_N);
  int nblocks = 1 + (ntot-1)/threadsPerBlock;
  int nthreads = threadsPerBlock;

  Update<<<nblocks, nthreads>>>(A, B, C);
  UpdateAtNode<<<nblocks, nthreads>>>(A, B, C);

  swap(sum, sum_cnt);
  swap(square, square_cnt);
  swap(P, P_cnt);
  swap(Theta, Theta_cnt);
  swap(Q00, Q00_cnt);
  swap(Q01, Q01_cnt);
  swap(Px, Px_cnt);
  swap(Py, Py_cnt);
  swap(area, area_cnt);
}
