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
    double save_zeta  = 0; swap(zeta,  save_zeta);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i)
      for(unsigned i=0; i<=npc; ++i) Update(i==0);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(zeta, save_zeta);
  }

  if(BC==5) ConfigureWalls(1);
  if(BC==6) ConfigureWalls(0);
}

void Model::Post()
{}

void Model::PreRunStats()
{
  // packing fraction
  {
    // total cell area
    double packing = 0.;
    for(unsigned n=0; n<nphases; ++n)
      packing += R*R;
    packing *= Pi;

    // divide by available area
    double area = 0;
    for(unsigned k=0; k<N; ++k)
      area += 1.-walls[k];
    packing /= area;

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
    if(abs(1.-area[n]/(Pi*R*R))>.2)
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
  const auto p  = phi[n][q];
  const auto dx = derivX(phi[n], sq);
  const auto dy = derivY(phi[n], sq);
  const auto a  = area[n];
  const auto a0 = Pi*R*R;
  // all-cells properties
  const auto ll = laplacian(phi[n], sq);
  const auto ls = laplacian(sum, s);

  // delta F / delta phi_i
  V[n][q] = (
      // CH term
      + gam*(8*p*(1-p)*(1-2*p)/lambda - 2*lambda*ll)
      // area conservation term
      - 4*mu/a0*(1-a/a0)*p
      // repulsion term
      + 4*kappa/lambda*p*(square[k]-p*p)
    );

  // isotropic part of the stress = sum_i delta F / delta phi_i
  const double pressure = (
      // CH term
      + gam*(8*(sum[k]-3*square[k]+2*thirdp[k])/lambda - 2*lambda*ls)
      // area conservation term
      - 4*mu/a0*(sum[k]-sumA[k]/a0)
      // repulsion term
      - 4*kappa/lambda*(sum[k]*square[k]-thirdp[k])
    );

  // passive force
  force_p[n][0] += -pressure*dx;
  force_p[n][1] += -pressure*dy;
  // contractility force
  force_c[n][0] += zeta*sumS00[k]*dx + zeta*sumS01[k]*dy;
  force_c[n][1] += zeta*sumS01[k]*dx - zeta*sumS00[k]*dy;
}

void Model::UpdateAtNode(unsigned n, unsigned q, bool store)
{
  const auto   k = GetIndexFromPatch(n, q);

  // compute dphi
  {
    const auto& sq = neighbors_patch[q];
    const auto dx  = derivX(phi[n], sq);
    const auto dy  = derivY(phi[n], sq);

    dphi[n][q] =
      // free energy
      - V[n][q]
      // advection term
      - velocity[n][0]*dx - velocity[n][1]*dy;
      ;
  }

  // store values
  if(store)
  {
    dphi_old[n][q] = dphi[n][q];
    phi_old[n][q]  = phi[n][q];
  }

  // predictor-corrector
  {
    double p = phi_old[n][q]
               + time_step*.5*(dphi[n][q] + dphi_old[n][q]);

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
  ReinitSumsAtNode(k);
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
  const coord com_grd { unsigned(round(com[n][0])), unsigned(round(com[n][1])) };
  const coord new_min = ( com_grd + Size - patch_margin ) % Size;
  const coord new_max = ( com_grd + patch_margin - coord {1u, 1u} ) % Size;
  coord displacement  = ( Size + new_min - patch_min[n] ) % Size;

  // I guess there is somehthing better than this...
  if(displacement[0]==Size[0]-1u) displacement[0] = patch_size[0]-1u;
  if(displacement[1]==Size[1]-1u) displacement[1] = patch_size[1]-1u;

  // update offset and patch location
  offset[n]    = ( offset[n] + patch_size - displacement ) % patch_size;
  patch_min[n] = new_min;
  patch_max[n] = new_max;
}

void Model::UpdateStructureTensorAtNode(unsigned n, unsigned q)
{
  // to be reintroduced correctly
  const auto& sq = neighbors_patch[q];
  const auto  dx = derivX(phi[n], sq);
  const auto  dy = derivY(phi[n], sq);

  S00[n] += -0.5*(dx*dx-dy*dy);
  S01[n] += -dx*dy;
}

void Model::UpdateSumsAtNode(unsigned n, unsigned q)
{
  const auto p = phi[n][q];
  const auto k = GetIndexFromPatch(n, q);

  sum_cnt[k]    += p;
  square_cnt[k] += p*p;
  thirdp_cnt[k] += p*p*p;
  sumA_cnt[k]   += p*area[n];
  sumS00_cnt[k] += p*S00[n];
  sumS01_cnt[k] += p*S01[n];
}

inline void Model::ReinitSumsAtNode(unsigned k)
{
  sum[k]    = 0;
  square[k] = 0;
  thirdp[k] = 0;
  sumA[k]   = 0;
  sumS00[k] = 0;
  sumS01[k] = 0;
}

#ifndef _CUDA_ENABLED

void Model::Update(bool store, unsigned nstart)
{
  // Compute induced force and passive velocity
  //
  // We need to loop over all nodes once before updating the phase fields
  // because the passive component of the velocity requires a integral
  // involving a derivative. This means that the field phi must be fully
  // updated first.

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    force_p[n] = {0., 0.};
    force_c[n] = {0., 0.};

    // update in restricted patch only
    for(unsigned q=0; q<patch_N; ++q)
      UpdateFieldsAtNode(n, q);

    // normalise and compute total forces and vel
    velocity[n]  = (force_p[n] + force_c[n])/xi;
  }

  // Predictor-corrector function for updating the phase fields
  //
  // The predictor corrector is such that it can be used many time in a row in
  // order to give you better precision, effectively giving higher order
  // approximations.

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    S00[n] = S01[n] = 0.;

    // only update fields in the restricted patch of field n
    for(unsigned q=0; q<patch_N; ++q)
    {
      UpdateAtNode(n, q, store);
      UpdateStructureTensorAtNode(n, q);
    }

    // update center of mass
    ComputeCoM(n);
    // update patch boundaries
    UpdatePatch(n);

    com_x[n] = com_y[n] = area[n] = 0.;
  }

  // Compute square and sum
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
      UpdateSumsAtNode(n, q);
  }

  // Reinit counters and swap to get correct values
  //
  // We use this construction because it is much faster with OpenMP: the single
  // threaded portion of the code consists only of these swaps!

  swap(area, area_cnt);
  swap(sum, sum_cnt);
  swap(square, square_cnt);
  swap(thirdp, thirdp_cnt);
  swap(sumA, sumA_cnt);
  swap(sumS00, sumS00_cnt);
  swap(sumS01, sumS01_cnt);
}

#endif
