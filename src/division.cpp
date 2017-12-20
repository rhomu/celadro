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

using namespace std;

void Model::Divide(unsigned n)
{
  // decrease counter
  division_counter[n] -= 1;

  if(division_counter[n]==nsubsteps*division_time)
  {
    cout << "START\n";
    R[n] *= division_growth;
  }

  if(division_counter[n]==0)
  {
    cout << "DIVIDE\n";

    // -------------------------------------------------------------------------

    // add the two daughters
    nphases += 2;

    // extend memory
    gam.resize(nphases, gam[n]);
    mu.resize(nphases, mu[n]);
    delta.resize(nphases, delta[n]);
    R.resize(nphases, R[n]/division_growth/2.);
    dR.resize(nphases, R[n]/2./nsubsteps/division_refract_time);
    phi.resize(nphases, vector<double>(patch_N, 0.));
    phi_old.resize(nphases, vector<double>(patch_N, 0.));
    V.resize(nphases, vector<double>(patch_N, 0.));
    potential.resize(nphases, vector<double>(patch_N, 0.));
    potential_old.resize(nphases, vector<double>(patch_N, 0.));
    division_counter.resize(nphases, 0.);
    area.resize(nphases, 0.);
    area_cnt.resize(nphases, 0.);
    patch_min.resize(nphases, {0, 0});
    patch_max.resize(nphases, Size);
    com.resize(nphases, {0., 0.});
    com_prev.resize(nphases, {0., 0.});
    pol.resize(nphases, {0., 0.});
    velp.resize(nphases, {0., 0.});
    velc.resize(nphases, {0., 0.});
    velf.resize(nphases, {0., 0.});
    com_x.resize(nphases, 0.);
    com_y.resize(nphases, 0.);
    c.resize(nphases, 0);
    S00.resize(nphases, 0.);
    S01.resize(nphases, 0.);
    S_order.resize(nphases, 0.);
    S_angle.resize(nphases, 0.);
    theta.resize(nphases, 0.);
    offset.resize(nphases, {0u, 0u});

    // create the cells
    const double theta0 = theta[n];
    const double length = .5*sqrt(area[n]/Pi);
    const vec<double, 2> d = {cos(theta0), sin(theta0)};
    AddCell(nphases-2, coord(com[n] + length*d));
    AddCell(nphases-1, coord(com[n] - length*d));

    // -------------------------------------------------------------------------
    // Relax daughter cells

    double save_c0      = 0.; swap(c0, save_c0);
    double save_zeta    = 0.; swap(zeta, save_zeta);
    double save_beta    = 0.; swap(beta, save_beta);
    double save_alpha   = 0.; swap(alpha, save_alpha);
    double save_S = 0.; swap(S, save_S);
    double save_J = 0.; swap(J, save_J);
    double save_D = 0.; swap(D, save_D);

    for(unsigned i=0; i<division_relax_time*nsubsteps; ++i)
    {
      for(unsigned i=0; i<=npc; ++i)
      {
        // add inverse potential to sum and square sum of phase fields
        PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
        for(unsigned q=0; q<patch_N; ++q)
        {
          const auto k = GetIndexFromPatch(n, q);
          const auto a = 1-phi[n][q];
          sum[k]      += a;
          square[k]   += a*a;
        }

        // only update last two cells only
        Update(i==0, nphases-2);
      }
    }

    swap(c0, save_c0);
    swap(zeta, save_zeta);
    swap(beta, save_beta);
    swap(alpha, save_alpha);
    swap(S, save_S);
    swap(J, save_J);
    swap(D, save_D);

    cout << "END\n";

    // -------------------------------------------------------------------------
    // destroy the old fat mama

    // swap last and n
    SwapCells(nphases-1, n);

    // shrink memory
    gam.pop_back();
    mu.pop_back();
    delta.pop_back();
    R.pop_back();
    phi.pop_back();
    phi_old.pop_back();
    V.pop_back();
    potential.pop_back();
    potential_old.pop_back();
    division_counter.pop_back();
    area.pop_back();
    area_cnt.pop_back();
    patch_min.pop_back();
    patch_max.pop_back();
    com.pop_back();
    com_prev.pop_back();
    pol.pop_back();
    velp.pop_back();
    velc.pop_back();
    velf.pop_back();
    com_x.pop_back();
    com_y.pop_back();
    c.pop_back();
    S00.pop_back();
    S01.pop_back();
    S_order.pop_back();
    S_angle.pop_back();
    theta.pop_back();
    offset.pop_back();

    nphases -= 1;
  }

  if(division_counter[n]<0)
  {
    R[n] += dR[n];

    if(unsigned(-division_counter[n]) == nsubsteps*division_refract_time)
      ResetDivisionCounter(n);
  }
  cout << n << ' ' << R[n] << endl;
}

void Model::ResetDivisionCounter(unsigned n)
{
  division_counter[n] = nsubsteps*(division_time + 1
      + random_exponential(division_rate));
}
