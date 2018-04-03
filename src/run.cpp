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
//    double save_alpha = 0; swap(alpha, save_alpha);
    double save_Dnem  = 0; swap(Dnem, save_Dnem);
    double save_Dpol  = 0; swap(Dpol, save_Dpol);
    double save_Jnem  = 0; swap(Jnem, save_Jnem);
    double save_Jpol  = 0; swap(Jpol, save_Jpol);
    double save_Kpol  = 0; swap(Kpol, save_Kpol);
    double save_Knem  = 0; swap(Knem, save_Knem);
    double save_Wnem  = 0; swap(Wnem, save_Wnem);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i)
      for(unsigned i=0; i<=npc; ++i) Update(i==0);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(zeta, save_zeta);
//    swap(alpha, save_alpha);
    swap(Jnem, save_Jnem);
    swap(Jpol, save_Jpol);
    swap(Dnem, save_Dnem);
    swap(Dpol, save_Dpol);
    swap(Kpol, save_Kpol);
    swap(Knem, save_Knem);
    swap(Wnem, save_Wnem);
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
      packing += R[n]*R[n];
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
  const auto  dxs  = derivX(sum, s);
  const auto  dys  = derivY(sum, s);
  const auto  dxp0 = derivX(P0, s);
  const auto  dyp0 = derivY(P0, s);
  const auto  dxp1 = derivX(P1, s);
  const auto  dyp1 = derivY(P1, s);
  // walls properties
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
  V[n][q]     = force;
  // passive force
  //force_p[n][0] += dx*force;
  //force_p[n][1] += dy*force;
  force_p[n][0] += C1*kappa*square[k]*dx;
  force_p[n][1] += C1*kappa*square[k]*dy;
  // contractility force
  force_c[n][0] += zeta*sumQ00[k]*dx + zeta*sumQ01[k]*dy;
  force_c[n][1] += zeta*sumQ01[k]*dx - zeta*sumQ00[k]*dy;
  // friction force
  force_f[n][0] += + f*alpha[n]*pol[n][0]*(dx*(dxs-dx)+dy*(dys-dy))
                   - f*alpha[n]*(dx*(dxp0-pol[n][0]*dx)+dy*(dyp0-pol[n][0]*dy))
                   - dyw*f_walls*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw);
  force_f[n][1] += + f*alpha[n]*pol[n][1]*(dx*(dxs-dx)+dy*(dys-dy))
                   - f*alpha[n]*(dx*(dxp1-pol[n][1]*dx)+dy*(dyp1-pol[n][1]*dy))
                   + dxw*f_walls*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw);
  // differential torques...
  {
    // ... overlap between cells
    const double ovlap = -(dx*(dxs-dx)+dy*(dys-dy));
    // ... polarisation torque
    const vec<double, 2> P = {P0[k]-phi[n][k]*pol[n][0], P1[k]-phi[n][k]*pol[n][1]};
    delta_theta_pol[n] += ovlap*atan2(P[0]*pol[n][1]-P[1]*pol[n][0],
                                       P[0]*pol[n][0]+P[1]*pol[n][1]);
    // ... nematic torque
    //const vec<double, 2> Q = {sumQ00[k]-phi[n][k]*Q00[n], sumQ01[k]-phi[n][k]*Q01[n]};
    //delta_theta_nem[n] += ovlap*atan2(Q[0]*Q01[n]-Q[1]*Q00[n],
    //                                  Q[0]*Q00[n]+Q[1]*Q01[n]);
  }
  // nematic torques
  tau[n]       += sumQ00[k]*Q01[n]-sumQ01[k]*Q00[n];
  vorticity[n] += U0[k]*dy-U1[k]*dx;
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
    // distance from the center of mass
    const auto dc = diff(vec<double, 2>(GetPosition(k)), com[n]);
    const auto el = pol[n].sq()<1e-6 ? 0. :
      - 2./(Pi*r*r)*delta[n]*p*(dc - pol[n]*(pol[n]*dc)/pol[n].sq()).sq();

    potential[n][q] = (
      // free energy term
      -.5*V[n][q]
      -.5*(
        + C1*gam[n]*p*(1.-p)*(1.-2.*p)
        - 2.*mu[n]*(1.-a/(Pi*r*r))*2.*p
        - 2.*gam[n]*ll
      )
      // elongation potential
      + el
      // advection term
      - velocity[n][0]*dx - velocity[n][1]*dy
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
  // tot force
  const double fn  = force_tot[n].abs();
  // force tensor
  const double F00 =   force_tot[n][0]*force_tot[n][0]
                     - force_tot[n][1]*force_tot[n][1];
  const double F01 = 2*force_tot[n][0]*force_tot[n][1];

  // update nematics and polarity
  if(store)
  {
    // euler-marijuana update
    theta_nem_old[n] = theta_nem[n] + sqrt_time_step*Dnem*random_normal();
    theta_pol_old[n] = theta_pol[n] + sqrt_time_step*Dpol*random_normal();
  }

  // nematics
  theta_nem[n] = theta_nem_old[n] - time_step*(
      + Knem*tau[n]
      + Jnem*fn*atan2(F00*Q01[n]-F01*Q00[n], F00*Q00[n]+F01*Q01[n]))
      + Wnem*vorticity[n];
  Q00[n] = Snem*cos(2*theta_nem[n]);
  Q01[n] = Snem*sin(2*theta_nem[n]);

  // polarisation
  theta_pol[n] = theta_pol_old[n] - time_step*(
      + Kpol*delta_theta_pol[n]
      + Jpol*fn*atan2(force_tot[n][0]*pol[n][1]-force_tot[n][1]*pol[n][0], force_tot[n]*pol[n]));
  pol[n] = { Spol*cos(theta_pol[n]), Spol*sin(theta_pol[n]) };
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
  coord displacement  = ( Size + new_min - patch_min[n]) % Size;

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

void Model::SquareAndSumAtNode(unsigned n, unsigned q)
{
  const auto p = phi[n][q];
  const auto k = GetIndexFromPatch(n, q);

  // we swap counters and values afterwards
  sum_cnt[k]    += p;
  square_cnt[k] += p*p;
  sumQ00_cnt[k] += p*Q00[n];
  sumQ01_cnt[k] += p*Q01[n];
  P0_cnt[k]     += p*pol[n][0];
  P1_cnt[k]     += p*pol[n][1];
  U0_cnt[k]     += p*velocity[n][0];
  U1_cnt[k]     += p*velocity[n][1];
}

inline void Model::ReinitSquareAndSumAtNode(unsigned k)
{
  sum[k]    = 0;
  square[k] = 0;
  sumQ00[k] = 0;
  sumQ01[k] = 0;
  P0[k]     = 0;
  P1[k]     = 0;
  U0[k]     = 0;
  U1[k]     = 0;
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
    force_p[n] = {0., 0.};
    force_c[n] = {0., 0.};
    force_f[n] = {0., 0.};
    delta_theta_pol[n] = 0;
    tau[n] = 0;
    vorticity[n] = 0;

    // update in restricted patch only
    for(unsigned q=0; q<patch_N; ++q)
      UpdateFieldsAtNode(n, q);

    // normalise and compute total forces and vel
    tau[n]      /= lambda;
    force_tot[n] = force_p[n] + force_c[n] + force_f[n];
    velocity[n]  = (force_tot[n] + alpha[n]*pol[n])/xi;
  }

  // 2) Predictor-corrector function for updating the phase fields
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
    // Update Q-tensor and polarisation
    UpdatePolarization(n, store);
    // update center of mass
    ComputeCoM(n);
    // update patch boundaries
    UpdatePatch(n);

    com_x[n] = com_y[n] = area[n] = 0.;
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
  swap(P0, P0_cnt);
  swap(P1, P1_cnt);
  swap(U0, U0_cnt);
  swap(U1, U1_cnt);
}

#endif
