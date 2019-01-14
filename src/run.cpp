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
    vector<double> save_alpha(nphases, 0); swap(alpha,  save_alpha);
    vector<double> save_zetaS(nphases, 0); swap(zetaS,  save_zetaS);
    vector<double> save_zetaQ(nphases, 0); swap(zetaQ,  save_zetaQ);
    vector<double> save_Dnem(nphases, 0); swap(Dnem,  save_Dnem);
    vector<double> save_Dpol(nphases, 0); swap(Dpol,  save_Dpol);
    vector<double> save_Jnem(nphases, 0); swap(Jnem,  save_Jnem);
    vector<double> save_Jpol(nphases, 0); swap(Jpol,  save_Jpol);
    vector<double> save_Knem(nphases, 0); swap(Knem,  save_Knem);
    vector<double> save_Wnem(nphases, 0); swap(Wnem,  save_Wnem);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i)
      for(unsigned i=0; i<=npc; ++i) Update(i==0);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(alpha, save_alpha);
    swap(zetaS, save_zetaS);
    swap(zetaQ, save_zetaQ);
    swap(Jnem, save_Jnem);
    swap(Jpol, save_Jpol);
    swap(Dnem, save_Dnem);
    swap(Dpol, save_Dpol);
    swap(Knem, save_Knem);
    swap(Wnem, save_Wnem);
  }

  if(BC==5 || BC==7) ConfigureWalls(1);
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

void Model::UpdateSumsAtNode(unsigned n, unsigned q)
{
  const auto k = GetIndexFromPatch(n, q);
  const auto p = phi[n][q];

  sum[k]     += p;
  square[k]  += p*p;
  thirdp[k]  += p*p*p;
  fourthp[k] += p*p*p*p;
  sumA[k]    += p*p*area[n];
  sumS00[k]  += p*S00[n]*zetaS[n];
  sumS01[k]  += p*S01[n]*zetaS[n];
  sumQ00[k]  += p*Q00[n]*zetaQ[n];
  sumQ01[k]  += p*Q01[n]*zetaQ[n];
  P0[k]      += p*polarization[n][0];
  P1[k]      += p*polarization[n][1];
  U0[k]      += p*velocity[n][0];
  U1[k]      += p*velocity[n][1];
}

void Model::UpdatePotAtNode(unsigned n, unsigned q)
{
  const auto  k  = GetIndexFromPatch(n, q);
  const auto& s  = neighbors[k];
  const auto& sq = neighbors_patch[q];

  const auto p  = phi[n][q];
  const auto a  = area[n];
  const auto ll = laplacian(phi[n], sq);
  const auto ls = laplacian(sum, s);

  const double internal = (
      // CH term
      + gam[n]*(8*p*(1-p)*(1-2*p)/lambda - 2*lambda*ll)
      // area conservation term
      - 4*mu[n]/a0[n]*(1-a/a0[n])*p
    );

  const double interactions = (
      // repulsion term
      + 2*kappa/lambda*p*(square[k]-p*p)
      // adhesion term
      - 2*omega*lambda*(ls-ll)
      // repulsion with walls
      + 2*wall_kappa/lambda*p*walls[k]*walls[k]
      // adhesion with walls
      - 2*wall_omega*lambda*walls_laplace[k]
    );

  // delta F / delta phi_i
  V[n][q] = internal + interactions;

  // pressure
  pressure[k] += p*interactions;
}

void Model::UpdateForcesAtNode(unsigned n, unsigned q)
{
  const auto  k  = GetIndexFromPatch(n, q);
  const auto& sq = neighbors_patch[q];

  const auto dx  = derivX(phi[n], sq);
  const auto dy  = derivY(phi[n], sq);

  stress_xx[k] = - pressure[k] - sumS00[k] - sumQ00[k];
  stress_yy[k] = - pressure[k] + sumS00[k] + sumQ00[k];
  stress_xy[k] = - sumS01[k] - sumQ01[k];

  Fpressure[n] += { pressure[k]*dx, pressure[k]*dy };
  Fshape[n]    += { sumS00[k]*dx + sumS01[k]*dy,
                    sumS01[k]*dx - sumS00[k]*dy };
  Fnem[n]      += { sumQ00[k]*dx + sumQ01[k]*dy,
                    sumQ01[k]*dx - sumQ00[k]*dy };

  // store derivatives
  phi_dx[n][q] = dx;
  phi_dy[n][q] = dy;

  // nematic torques
  tau[n]       += sumQ00[k]*Q01[n] - sumQ01[k]*Q00[n];
  vorticity[n] += U0[k]*dy - U1[k]*dx;
}

void Model::UpdatePhaseFieldAtNode(unsigned n, unsigned q, bool store)
{
  const auto k = GetIndexFromPatch(n, q);

  // compute dphi
  dphi[n][q] =
    // free energy
    - V[n][q]
    // advection term
    - velocity[n][0]*phi_dx[n][q] - velocity[n][1]*phi_dy[n][q];
    ;

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

    com_x[n] += com_x_table[GetXPosition(k)]*p;
    com_y[n] += com_y_table[GetYPosition(k)]*p;
    area[n]  += p*p;
  }

  // reinit values: we do reinit values here for the simple reason that it is
  // faster than having a supplementary loop afterwards. There is a race
  // condition in principle here but since we are setting evth back to 0 it
  // should be fine. Note that this should be done before the patches are updated
  ReinitSumsAtNode(k);
}

void Model::UpdateNematic(unsigned n, bool store)
{
  // euler-marijuana update
  if(store)
    theta_nem_old[n] = theta_nem[n] + sqrt_time_step*Dnem[n]*random_normal();

  double F00 = 0, F01 = 0;
  switch(align_nematic_to)
  {
    case 0:
    {
      const auto ff = velocity[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
    }
    case 1:
    {
      const auto ff = Fpressure[n];
      F00 =   ff[0]*ff[0]
            - ff[1]*ff[1];
      F01 = 2*ff[0]*ff[1];
      break;
    }
    case 2:
      F00 = S00[n];
      F01 = S01[n];
      break;
  }
  const auto strength = pow(F01*F01 + F00*F00, 0.25);

  theta_nem[n] = theta_nem_old[n] - time_step*(
      + Knem[n]*tau[n]
      + Jnem[n]*strength*atan2(F00*Q01[n]-F01*Q00[n], F00*Q00[n]+F01*Q01[n]))
      + Wnem[n]*vorticity[n];
  Q00[n] = Snem[n]*cos(2*theta_nem[n]);
  Q01[n] = Snem[n]*sin(2*theta_nem[n]);
}

void Model::UpdatePolarization(unsigned n, bool store)
{
  // euler-marijuana update
  if(store)
    theta_pol_old[n] = theta_pol[n] + sqrt_time_step*Dpol[n]*random_normal();

  vec<double, 2> ff = {0, 0};
  switch(align_polarization_to)
  {
    case 0:
      ff = velocity[n];
      break;
    case 1:
      ff = Fpressure[n];
      break;
  }

  theta_pol[n] = theta_pol_old[n] - time_step*(
      + Jpol[n]*ff.abs()*atan2(ff[0]*polarization[n][1]-ff[1]*polarization[n][0], ff*polarization[n]));
  polarization[n] = { Spol[n]*cos(theta_pol[n]), Spol[n]*sin(theta_pol[n]) };
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
  const auto  dx = phi_dx[n][q];
  const auto  dy = phi_dy[n][q];

  S00[n] += -0.5*(dx*dx-dy*dy);
  S01[n] += -dx*dy;
}

void Model::ReinitSumsAtNode(unsigned k)
{
  sum[k] = 0;
  square[k] = 0;
  thirdp[k] = 0;
  fourthp[k] = 0;
  sumA[k] = 0;
  sumS00[k] = 0;
  sumS01[k] = 0;
  sumQ00[k] = 0;
  sumQ01[k] = 0;
  pressure[k] = 0;
  U0[k] = 0;
  U1[k] = 0;
}

#ifndef _CUDA_ENABLED

void Model::Update(bool store, unsigned nstart)
{
  // Compute all global sums
  for(unsigned n=nstart; n<nphases; ++n)
  {
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      UpdateSumsAtNode(n, q);
  }

  // Compute stresses
  //
  // We need another loop because the pressure involves a double sum over all
  // the cells.
  for(unsigned n=nstart; n<nphases; ++n)
  {
    // update only patch (in parallel, each node to a different core)
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
      UpdatePotAtNode(n, q);
  }

  // Compute induced force and passive velocity
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    Fpol[n] = Fshape[n] = Fpressure[n] = {0, 0};
    tau[n] = vorticity[n] = 0;

    // update in restricted patch only
    for(unsigned q=0; q<patch_N; ++q)
      UpdateForcesAtNode(n, q);

    // normalise and compute total forces and vel
    tau[n]     /= lambda;
    Fpol[n]     = alpha[n]*polarization[n];
    velocity[n] = (Fpressure[n] + Fnem[n] + Fshape[n] + Fpol[n])/xi[n];
  }

  // Predictor-corrector function for updating the phase fields
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=nstart; n<nphases; ++n)
  {
    com_x[n] = com_y[n] = S00[n] = S01[n] = area[n] = 0;

    // only update fields in the restricted patch of field n
    for(unsigned q=0; q<patch_N; ++q)
    {
      UpdatePhaseFieldAtNode(n, q, store);
      UpdateStructureTensorAtNode(n, q);
    }

    // update polarisation
    UpdatePolarization(n, store);
    // update Q-tensor
    UpdateNematic(n, store);
    // update center of mass
    ComputeCoM(n);
    // update patch boundaries
    UpdatePatch(n);
  }
}

#endif
