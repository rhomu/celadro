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
#include "random.hpp"
#include "derivatives.hpp"
#include "tools.hpp"

using namespace std;

/** Needa store? (yes this is against my will) */
bool store;

void Model::Pre()
{
  // we make the system relax (without activity)
  if(relax_time>0)
  {
    double save_c0 = 0.; swap(c0, save_c0);
    double save_zeta = 0.; swap(zeta, save_zeta);
    double save_beta = 0.; swap(beta, save_beta);
    double save_alpha = 0.; swap(alpha, save_alpha);
    double save_division_rate = 0.; swap(save_division_rate, division_rate);

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    for(unsigned i=0; i<relax_time*nsubsteps; ++i) Step();

    if(relax_nsubsteps) swap(nsubsteps, relax_nsubsteps);

    swap(c0, save_c0);
    swap(zeta, save_zeta);
    swap(beta, save_beta);
    swap(alpha, save_alpha);
    swap(save_division_rate, division_rate);
  }
}

void Model::Post()
{}

void Model::PreRunStats()
{
  // packing fraction
  {
    double packing = nphases*Pi*R*R;
    if(BC<=2)
      packing/= (birth_bdries[1]-birth_bdries[0])
               *(birth_bdries[3]-birth_bdries[2]);
    if(BC==3) packing /= Pi*Size/4.;

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
  for(const auto a : area)
    if(abs(1.-a/C2)>.2)
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

void Model::UpdateFieldsAtNode(unsigned n, unsigned k)
{
  const auto& s = neighbors[k];

  // cell properties
  const auto& p  = phi[n][k];
  const auto  l  = laplacian(phi[n], s);
  const auto  dx = derivX(phi[n], s);
  const auto  dy = derivY(phi[n], s);
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
        + omega*(ls-l)
        + wall_omega*lw
        )
      );

  // potential
  V[n][k]     = force;
  // passive force
  velp[n][0] += dx*force;
  velp[n][1] += dy*force;
  // contractility force
  velc[n][0] += ( (P[k]+zeta*Q00[k])*dx + zeta*Q01[k]*dy );
  velc[n][1] += ( zeta*Q01[k]*dx + (P[k]-zeta*Q00[k])*dy );
  // friction force
  velf[n][0] += + f*alpha/xi*pol[n][0]*(dx*(dxs-dx)+dy*(dys-dy))
                - f*alpha/xi*(dx*(dxp0-pol[n][0]*dx)+dy*(dyp0-pol[n][0]*dy))
                //+ f_walls*alpha/xi*pol[n][0]*(dx*dxw+dy*dyw);
                - dyw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw);
  velf[n][1] += + f*alpha/xi*pol[n][1]*(dx*(dxs-dx)+dy*(dys-dy))
                - f*alpha/xi*(dx*(dxp1-pol[n][1]*dx)+dy*(dyp1-pol[n][1]*dy))
                //+ f_walls*alpha/xi*pol[n][1]*(dx*dxw+dy*dyw);
                + dxw*f_walls/xi*(pol[n][1]*dxw-pol[n][0]*dyw)*(dx*dxw+dy*dyw);

  // these are different alignment torques and must be cleaned up once we decide
  // which one is the best
  //
  // alignment torque (orientational)
  //const auto dxQ00 = derivX(Q00, s);
  //const auto dxQ01 = derivX(Q01, s);
  //const auto dyQ00 = derivY(Q00, s);
  //const auto dyQ01 = derivY(Q01, s);
  //torque[n] += - 4*J1*(Q00[n]*(dx*dxQ01+dy*dyQ01)-Q01[n]*(dx*dxQ00+dy*dyQ00));
  //
  // alignment torque (directional)
  //const auto  dxt  = derivX(Theta, s);
  //const auto  dyt  = derivY(Theta, s);
  //torque[n] -= J1*(dx*dxt+dy*dyt - theta[n]*(dx*dxs+dy*dys));
  //
  // alignment torque (shear stress)
  //torque[n] -= J1*zeta*(2*Q00[n]*dx*dy + Q01[n]*(dx*dx-dy*dy));
}

void Model::UpdateAtNode(unsigned n, unsigned k)
{
  // compute potential
  {
    const auto& s = neighbors[k];

    const auto  p  = phi[n][k];
    const auto  a  = area[n];
    const auto  dx = derivX(phi[n], s);
    const auto  dz = derivY(phi[n], s);
    const auto  l  = laplacian(phi[n], s);

    potential[n][k] = (
      // free energy term
      -.5*V[n][k]
      -.5*(
        + C1*gam[n]*p*(1.-p)*(1.-2.*p)
        - 2.*mu[n]*(1.-a/C2)*2.*p
        - 2.*gam[n]*l
      )
      // advection term
      - (alpha*pol[n][0]+velp[n][0]+velc[n][0]+velf[n][0])*dx/xi
      - (alpha*pol[n][1]+velp[n][1]+velc[n][1]+velf[n][1])*dz/xi
      );
  }

  // store values
  if(store)
  {
    potential_old[n][k] = potential[n][k];
    phi_old[n][k]       = phi[n][k];
  }

  // predictor-corrector
  {
    double p = phi_old[n][k]
               + time_step*.5*(potential[n][k] + potential_old[n][k]);

    // update for next call
    phi[n][k]    = p;
    com_x[n]    += com_x_table[GetXPosition(k)]*p;
    com_y[n]    += com_y_table[GetYPosition(k)]*p;
    area_cnt[n] += p*p;
  }
}

void Model::UpdatePolarization(unsigned n)
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
  pol[n] = {cos(theta[n]), sin(theta[n])};

  // dynamics of the contractility: needs some cleaning up once settled
  //
  //c[n]  += -time_step*(c[n]-c0)/tauc + beta*c0*(area_cnt[n]-area[n]);
  //c[n]  += time_step*( -(c[n]-c0) - beta*(1.-area_cnt[n]/C2/.88171))/tauc;
  //c[n]  -= time_step*(c[n]-beta*c0*area_cnt[n]/C2)/tauc;
}

void Model::ComputeCoM(unsigned n)
{
  // the strategy to deal with the periodic boundary conditions is to compute
  // all the integrals in Fourier space and come back at the end. This way the
  // periodicity of the domain is automatically taken into account.
  const auto mx = arg(com_x[n]/static_cast<double>(Size)) + Pi;
  const auto my = arg(com_y[n]/static_cast<double>(Size)) + Pi;
  com[n] = { mx/2./Pi*LX, my/2./Pi*LY };
}

void Model::UpdateWindow(unsigned n)
{
  // update walls
  domain_min[n][0] = (static_cast<unsigned>(round(com[n][0])) + LX - margin)%LX;
  domain_min[n][1] = (static_cast<unsigned>(round(com[n][1])) + LY - margin)%LY;
  domain_max[n][0] = (static_cast<unsigned>(round(com[n][0])) + margin)%LX;
  domain_max[n][1] = (static_cast<unsigned>(round(com[n][1])) + margin)%LY;
}

void Model::UpdateStructureTensorAtNode(unsigned n, unsigned k)
{
  const auto& s = neighbors[k];

  const auto  p  = phi[n][k];
  const auto  dx = 2*p*derivX(phi[n], s);
  const auto  dy = 2*p*derivY(phi[n], s);
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

void Model::SquareAndSumAtNode(unsigned n, unsigned k)
{
  const auto p = phi[n][k];

  // we swap counters and values afterwards
  sum_cnt[k]    += p;
  square_cnt[k] += p*p;
  P_cnt[k]      -= p*c[n];
  Q00_cnt[k]    -= p*.5*(pol[n][0]*pol[n][0]-pol[n][1]*pol[n][1]);
  Q01_cnt[k]    -= p*pol[n][0]*pol[n][1];
  Px_cnt[k]     += p*pol[n][0];
  Py_cnt[k]     += p*pol[n][1];
  Theta_cnt[k]  += p*theta[n];

  // coupling to shape
  //Q00_cnt[k] += p*p*(-c[n]+zeta*S_order[n]*cos(2*S_angle[n]));
  //Q01_cnt[k] += p*p*(     +zeta*S_order[n]*sin(2*S_angle[n]));
  //Q00_cnt[k] += p*p*(-c[n]-zeta*S_order[n]*cos(2*S_angle[n]));
}

inline void Model::ReinitSquareAndSumAtNode(unsigned k)
{
  // reinit values
  sum[k]    = 0;
  square[k] = 0;
  P[k]      = 0;
  Theta[k]  = 0;
  Q00[k]    = 0;
  Q01[k]    = 0;
  Px[k]     = 0;
  Py[k]     = 0;
}

/*
void Divide(unsigned n)
{
  // work in progress

  const double r = random_real();

  if(r<nsubsteps*division_rate)
  //if(not pine--)
  {
    cout << "DIVIDE" << endl;
    // add a new cell and extend memory
    const auto m = nphases++;
    InitializeFields();
    gam.resize(nphases, gam[n]);
    mu.resize(nphases, mu[n]);
    theta[m] = theta[n];
    // pick random axis
    //const double theta0 = 2*Pi*r/nsubsteps/division_rate;
    // set division axis to be perp to the polarization
    const double theta0 = theta[n] + Pi;

    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned k=0; k<Size; ++k)
    {
      const auto x = wrap(diff(unsigned(com[n][0]), GetXPosition(k)), LX);
      const auto y = wrap(diff(unsigned(com[n][1]), GetXPosition(k)), LY);
      const auto theta = atan2(y, x);

      // distance from cutting line
      const auto d = sqrt(x*x+y*y)*sin(theta-theta0);
      // current size of the cell (more or less)
      const auto l = 30;//5*sqrt(area[n]/Pi);

      const auto p = phi[n][k];
      phi[n][k] = .5*(1.+tanh( d/l))*p;
      phi[m][k] = .5*(1.-tanh(-d/l))*p;
    }
  }
}*/

void Model::Update()
{
  // 1) Compute induced force and passive velocity
  //
  // We need to loop over all nodes once before updating the phase fields
  // because the passive component of the velocity requires a integral
  // involving a derivative. This means that the field phi must be fully
  // updated first.

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=0; n<nphases; ++n)
  {
    velp[n] = {0., 0.};
    velc[n] = {0., 0.};
    velf[n] = {0., 0.};
    //torque[n] = 0.;

    // update in restricted domain only
    UpdateDomain(&Model::UpdateFieldsAtNode, n);
  }

  // 2) Predictor-corrector function for updating the phase fields
  //
  // The predictor corrector is such that it can be used many time in a row in
  // order to give you better precision, effectively giving higher order
  // approximations.

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=0; n<nphases; ++n)
  {
    // only update fields in the restricted domain of field n
    UpdateDomain(&Model::UpdateAtNode, n);
    // because the polarisation dynamics is first
    // order (euler-maruyama) we need to update only
    // once per predictor-corrector step
    if(store) UpdatePolarization(n);
    // update center of mass
    ComputeCoM(n);
    // update domain walls
    if(tracking) UpdateWindow(n);
    // update structure tensor
    UpdateDomain(&Model::UpdateStructureTensorAtNode, n);
    // and get shape
    ComputeShape(n);

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

  for(unsigned n=0; n<nphases; ++n)
    // update only domain (in parallel, each node to a different core)
    UpdateDomainP(&Model::SquareAndSumAtNode, n);

  // 4) Reinit counters and swap to get correct values
  //
  // We use this construction because it is much faster with OpenMP: the single
  // threaded portion of the code consists only of these swaps!

  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned k=0; k<Size; ++k) ReinitSquareAndSumAtNode(k);

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

void Model::Step()
{
  // first sweeps produces estimate of values
  store = true; // ok this is really bad :-(
  Update();
  store = false;

  // subsequent sweeps produce corrected values
  for(unsigned i=0; i<npc; ++i)
    Update();

  // division
  //const auto m = nphases; // this is needed because we might increment nphases
  //for(unsigned n=0; n<m; ++n) Divide(n);

  // compute center-of-mass velocity
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned n=0; n<nphases; ++n)
  {
    vel[n]      = { (com[n][0]-com_prev[n][0])/time_step,
                    (com[n][1]-com_prev[n][1])/time_step };
    com_prev[n] = com[n];
  }
}

// =============================================================================
// Helper functions

template<typename Ret, typename ...Args>
void Model::UpdateSubDomain(Ret (Model::*fun)(unsigned, unsigned, Args...),
                                unsigned n,
                                unsigned m0, unsigned m1,
                                unsigned M0, unsigned M1,
                                Args&&... args)
{
  // only update on the subregion
  for(unsigned i=m0; i<M0; ++i)
    for(unsigned j=m1; j<M1; ++j)
      // if you want to look it up, this is called a pointer
      // to member function and is an obscure C++ feature...
      (this->*fun)(n, GetDomainIndex(i, j),
                std::forward<Args>(args)...);
}

template<typename Ret, typename ...Args>
void Model::UpdateSubDomainP(Ret (Model::*fun)(unsigned, unsigned, Args...),
                                 unsigned n,
                                 unsigned m0, unsigned m1,
                                 unsigned M0, unsigned M1,
                                 Args&&... args)
{
  // same but with openmp
  PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
  for(unsigned k=0; k<(M0-m0)*(M1-m1); ++k)
      // if you want to look it up, this is called a pointer
      // to member function and is an obscure C++ feature...
      (this->*fun)(n, GetDomainIndex(m0+k%(M0-m0), m1+k/(M0-m0)),
                std::forward<Args>(args)...);
}

template<typename Ret, typename ...Args>
void Model::UpdateDomainP(Ret (Model::*fun)(unsigned, unsigned, Args...),
                              unsigned n, Args&&... args)
{
  if(domain_min[n][0]>=domain_max[n][0] and
     domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the corners
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], LX, LY);
    UpdateSubDomainP(fun, n, 0u, 0u, domain_max[n][0], domain_max[n][1]);
    UpdateSubDomainP(fun, n, domain_min[n][0], 0u, LX, domain_max[n][1]);
    UpdateSubDomainP(fun, n, 0u, domain_min[n][1], domain_max[n][0], LY);
  }
  else if(domain_min[n][0]>=domain_max[n][0])
  {
    // domain is across the left/right border
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], LX, domain_max[n][1]);
    UpdateSubDomainP(fun, n, 0u, domain_min[n][1], domain_max[n][0], domain_max[n][1]);
  }
  else if(domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the up/down border
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], LY);
    UpdateSubDomainP(fun, n, domain_min[n][0], 0u, domain_max[n][0], domain_max[n][1]);
  }
  else
    // domain is in the middle
    UpdateSubDomainP(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], domain_max[n][1]);
}

template<typename Ret, typename ...Args>
void Model::UpdateDomain(Ret (Model::*fun)(unsigned, unsigned, Args...),
                             unsigned n, Args&&... args)
{
  if(domain_min[n][0]>=domain_max[n][0] and
     domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the corners
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], LX, LY);
    UpdateSubDomain(fun, n, 0u, 0u, domain_max[n][0], domain_max[n][1]);
    UpdateSubDomain(fun, n, domain_min[n][0], 0u, LX, domain_max[n][1]);
    UpdateSubDomain(fun, n, 0u, domain_min[n][1], domain_max[n][0], LY);
  }
  else if(domain_min[n][0]>=domain_max[n][0])
  {
    // domain is across the left/right border
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], LX, domain_max[n][1]);
    UpdateSubDomain(fun, n, 0u, domain_min[n][1], domain_max[n][0], domain_max[n][1]);
  }
  else if(domain_min[n][1]>=domain_max[n][1])
  {
    // domain is across the up/down border
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], LY);
    UpdateSubDomain(fun, n, domain_min[n][0], 0u, domain_max[n][0], domain_max[n][1]);
  }
  else
    // domain is in the middle
    UpdateSubDomain(fun, n, domain_min[n][0], domain_min[n][1], domain_max[n][0], domain_max[n][1]);
}
