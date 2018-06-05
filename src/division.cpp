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
    cout << "START - phase " << n << "\n";
    R[n] *= division_growth;
  }

  if(division_counter[n]==0)
  {
    cout << "DIVIDE - phase " << n << "\n";

    // -------------------------------------------------------------------------

    // add the two daughters
    SetCellNumber(nphases + 2);

    // create the cells
    const double theta0 = atan2(S01[n], S00[n])/2.;
    const double length = .5*sqrt(area[n]/Pi);
    const auto d = length * vec<double, 2> {cos(theta0), sin(theta0)};
    AddCell(nphases-2, coord(com[n] + d));
    AddCell(nphases-1, coord(com[n] - d));

    // -------------------------------------------------------------------------
    // Relax daughter cells

    vector<double> save_alpha(nphases, 0); swap(alpha, save_alpha);
    double save_zeta  = 0; swap(zeta,  save_zeta);
    double save_Dnem  = 0; swap(Dnem,  save_Dnem);
    double save_Dpol  = 0; swap(Dpol,  save_Dpol);
    double save_Jnem  = 0; swap(Jnem,  save_Jnem);
    double save_Jpol  = 0; swap(Jpol,  save_Jpol);
    double save_Kpol  = 0; swap(Kpol,  save_Kpol);
    double save_Knem  = 0; swap(Knem,  save_Knem);
    double save_Wnem  = 0; swap(Wnem,  save_Wnem);

    // add inverse potential to sum and square sum of phase fields
    /*
    PRAGMA_OMP(omp parallel for num_threads(nthreads) if(nthreads))
    for(unsigned q=0; q<patch_N; ++q)
    {
      const auto k = GetIndexFromPatch(n, q);
      const auto a = 1-phi[n][q];
      sum[k]      += a;
      square[k]   += a*a;
    }
    */

    for(unsigned i=0; i<division_relax_time*nsubsteps; ++i)
      for(unsigned i=0; i<=npc; ++i)
        // only update last two cells only
        Update(i==0, nphases-2);

    swap(alpha, save_alpha);
    swap(zeta, save_zeta);
    swap(Jnem, save_Jnem);
    swap(Jpol, save_Jpol);
    swap(Dnem, save_Dnem);
    swap(Dpol, save_Dpol);
    swap(Kpol, save_Kpol);
    swap(Knem, save_Knem);
    swap(Wnem, save_Wnem);

    cout << "END - phase " << n << "\n";

    // -------------------------------------------------------------------------
    // destroy the old fat mama

    // swap last and n
    SwapCells(nphases-1, n);

    // shrink memory
    SetCellNumber(nphases-1);
  }

  if(division_counter[n]<0)
  {
    R[n] += dR[n];

    if(unsigned(-division_counter[n]) == nsubsteps*division_refract_time)
      ResetDivisionCounter(n);
  }
}

void Model::ResetDivisionCounter(unsigned n)
{
  division_counter[n] = nsubsteps*(division_time + 1
      + random_exponential(division_rate));
}
