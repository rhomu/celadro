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
#include "derivatives.hpp"
#include "tools.hpp"

using namespace std;

void Model::Pre()
{
  // we make the system relax (without activity)
  if(relax_time>0)
  {
    double save_zeta = 0.; swap(zeta, save_zeta);
    double save_beta = 0.; swap(beta, save_beta);
    double save_alpha = 0.; swap(alpha, save_alpha);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i)
      for(unsigned i=0; i<=npc; ++i) Update(i==0);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(zeta, save_zeta);
    swap(beta, save_beta);
    swap(alpha, save_alpha);
  }
}

void Model::Post()
{}

void Model::PreRunStats()
{
  // packing fraction
  {
    double packing = 0.;
    for(unsigned n=0; n<nphases; ++n)
      packing += R[n]*R[n];
    packing *= Pi;

    if(BC<=2)
      packing/= (birth_bdries[1]-birth_bdries[0])
               *(birth_bdries[3]-birth_bdries[2]);
    if(BC==3) packing /= Pi*N/4.;

    cout << "Packing fraction = " << packing << endl;
  }
}

void Model::RuntimeStats()
{
  // TBD
}

void Model::RuntimeChecks()
{
  // check that the area is more or less conserved (20%)
  for(unsigned n=0; n<nphases; ++n)
    if(abs(1.-area[n]/(Pi*R[n]*R[n]))>.2)
      throw warning_msg("area is not conserved.");

  for(unsigned n=0; n<nphases; ++n)
  {
    // check that the cells are not leaking, i.e. that at least 90% of the
    // contributions to the area comes from inside the cell (>1/2).
    double a = 0.;
    // compute area of points outside the cell (<1/2)
    for(const auto v : phi[n]) if(v<.5) a += v*v;
    // check that it is less than 5%
    if(a/area[n]>.90)
      throw warning_msg("your cells are leaking!");

    // check that the phase fields stay between 0 and 1
    for(const auto& p : phi[n])
      if(p<-0.5 or p>1.5)
        throw warning_msg("phase-field is not in [0,1]!");
  }
}

void Model::UpdateFieldsAtNode(unsigned n, unsigned q)
{
  const auto   k = GetIndexFromPatch(n, q);
  const auto&  s = neighbors[k];
  const auto& sq = neighbors_patch[q];

  // cell properties
  const auto& p  = phi[n][q];
  const auto  ll = laplacian(phi[n], sq);
  const auto  dx = derivX(phi[n], sq);
  const auto  dy = derivY(phi[n], sq);
  // all-cells properties
  const auto  ls   = laplacian(sum, s);
  // walls properties
  const auto  lw = walls_laplace[k];

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
  V[n][q]     = force;
  // passive force
  velp[n][0] += dx*force;
  velp[n][1] += dy*force;
  // contractility force
  velc[n][0] += ( (P+zeta*sumQ00[k])*dx + zeta*sumQ01[k]*dy );
  velc[n][1] += ( zeta*sumQ01[k]*dx + (P-zeta*sumQ00[k])*dy );
  // difference in contractility
  const double r = dx*dx+dy*dy;
  if(r>1e-6)
  {
    deltaQ00[n] += (sumQ00[k] - 2*phi[n][q]*Q00[n])*phi[n][q]*(Q00[n]*(dx*dx-dy*dy)+2*Q01[n]*dx*dy)/r;
    deltaQ01[n] += (sumQ01[k] - 2*phi[n][q]*Q01[n])*phi[n][q]*(Q00[n]*(dx*dx-dy*dy)+2*Q01[n]*dx*dy)/r;
  }
}

void Model::UpdateAtNode(unsigned n, unsigned q, bool store)
{
  const auto   k = GetIndexFromPatch(n, q);
  const auto& sq = neighbors_patch[q];

  // compute potential
  {
    const auto  p  = phi[n][q];
    const auto  a  = area[n];
    const auto  r  = R[n];
    const auto  dx = derivX(phi[n], sq);
    const auto  dy = derivY(phi[n], sq);
    const auto  ll = laplacian(phi[n], sq);

    potential[n][q] = (
      // free energy term
      -.5*V[n][q]
      -.5*(
        + C1*gam[n]*p*(1.-p)*(1.-2.*p)
        - 2.*mu[n]*(1.-a/(Pi*r*r))*2.*p
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
    potential_old[n][q] = potential[n][q];
    phi_old[n][q]       = phi[n][q];
  }

  // predictor-corrector
  {
    double p = phi_old[n][q]
               + time_step*.5*(potential[n][q] + potential_old[n][q]);

    // update for next call
    phi[n][q]    = p;
    com_x[n]    += com_x_table[GetXPosition(k)]*p;
    com_y[n]    += com_y_table[GetYPosition(k)]*p;
    area_cnt[n] += p*p;
  }

  // reinit values: we do reinit values here for the simple reason that it is
  // faster than having a supplementary loop afterwards. There is a race
  // condition in principle here but since we are setting evth back to 0 it
  // should be fine
  ReinitSquareAndSumAtNode(k);
}

void Model::UpdatePolarization(unsigned n, bool store)
{
  if(store)
  {
    Q00_old[n] = Q00[n] + sqrt(time_step)*D*random_normal();
    Q01_old[n] = Q01[n] + sqrt(time_step)*D*random_normal();
  }

  const double tr = 1 - Q00[n]*Q00[n] - Q01[n]*Q01[n];
  Q00[n] = Q00_old[n] + time_step*.5*(C*tr*Q00[n] + K*deltaQ00[n]);
  Q01[n] = Q01_old[n] + time_step*.5*(C*tr*Q01[n] + K*deltaQ01[n]);
}

void Model::ComputeCoM(unsigned n)
{
  // the strategy to deal with the periodic boundary conditions is to compute
  // all the integrals in Fourier space and come back at the end. This way the
  // periodicity of the domain is automatically taken into account.
  const auto mx = arg(com_x[n]/static_cast<double>(N)) + Pi;
  const auto my = arg(com_y[n]/static_cast<double>(N)) + Pi;
  com[n] = { mx/2./Pi*Size[0], my/2./Pi*Size[1] };
}

void Model::UpdatePatch(unsigned n)
{
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
}

void Model::UpdateStructureTensorAtNode(unsigned n, unsigned q)
{
  // to be reintroduced correctly
  const auto&   sq = neighbors_patch[q];

  const auto  p  = phi[n][q];
  const auto  dx = 2*p*derivX(phi[n], sq);
  const auto  dy = 2*p*derivY(phi[n], sq);

  S00[n] += 0.5*(dx*dx-dy*dy);
  S01[n] += dx*dy;
}

void Model::ComputeShape(unsigned n)
{
  // shape: we remap the 2x2 traceless symmetric matrix to polar coord for ease
  // of manipulation
  S_order[n] = sqrt(S00[n]*S00[n]+S01[n]*S01[n]);
  S_angle[n] = atan2(S01[n], S00[n])/2.;
}

void Model::SquareAndSumAtNode(unsigned n, unsigned q)
{
  const auto p = phi[n][q];
  const auto k = GetIndexFromPatch(n, q);

  // we swap counters and values afterwards
  sum_cnt[k]    += p;
  square_cnt[k] += p*p;
  sumQ00_cnt[k] += p*Q00[n];
  sumQ01_cnt[k] += p*Q01[n];
}

inline void Model::ReinitSquareAndSumAtNode(unsigned k)
{
  sum[k]    = 0;
  square[k] = 0;
  sumQ00[k] = 0;
  sumQ01[k] = 0;
}

#ifndef _CUDA_ENABLED

void Model::Update(bool store, unsigned nstart)
{
  // 1) Compute induced force and passive velocity
  //
  // We need to loop over all nodes once before updating the phase fields
  // because the passive component of the velocity requires a integral
  // involving a derivative. This means that the field phi must be fully
  // updated first.

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    velp[n] = {0., 0.};
    velc[n] = {0., 0.};
    velf[n] = {0., 0.};
    deltaQ00[n] = 0;
    deltaQ01[n] = 0;

    // update in restricted patch only
    for(unsigned q=0; q<patch_N; ++q)
      UpdateFieldsAtNode(n, q);
  }

  // 2) Predictor-corrector function for updating the phase fields
  //
  // The predictor corrector is such that it can be used many time in a row in
  // order to give you better precision, effectively giving higher order
  // approximations.

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    // only update fields in the restricted patch of field n
    for(unsigned q=0; q<patch_N; ++q)
      UpdateAtNode(n, q, store);
    // because the polarisation dynamics is first
    // order (euler-maruyama) we need to update only
    // once per predictor-corrector step
    /*if(store)*/ UpdatePolarization(n, store);
    // update center of mass
    ComputeCoM(n);
    // get shape
    ComputeShape(n);
    // update patch boundaries
    UpdatePatch(n);

    // reinit for next round
    com_x[n] = com_y[n] = area[n] = 0.;
    S00[n] = S01[n] = 0.;
  }

  // 3) Compute square and sum
  //
  // We need yet another loop here for parallelization, because we can not send
  // each phi to a different core when computing the square and the sum of all
  // phase fields. This is much faster than using an atomic portion in the
  // previous loop.

  for(unsigned n=nstart; n<nphases; ++n)
  {
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      SquareAndSumAtNode(n, q);
  }

  // 5) Reinit counters and swap to get correct values
  //
  // We use this construction because it is much faster with OpenMP: the single
  // threaded portion of the code consists only of these swaps!

  swap(area, area_cnt);
  swap(sum, sum_cnt);
  swap(square, square_cnt);
  swap(sumQ00, sumQ00_cnt);
  swap(sumQ01, sumQ01_cnt);
}

#endif
