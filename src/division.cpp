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

  if(division_counter[n]==nsubsteps*division_time[n])
  {
    cout << "START - phase " << n << "\n";
    R[n] *= division_growth;
  }

  if(division_counter[n]==0)
  {
    cout << "DIVIDE - phase " << n << "\n";

    // -------------------------------------------------------------------------

    // add the two daughters
    nphases += 2;

    gam.resize(nphases, gam[n]);
    mu.resize(nphases, mu[n]);
    R.resize(nphases, R[n]/division_growth);
    xi.resize(nphases, xi[n]);
    alpha.resize(nphases, alpha[n]);
    beta.resize(nphases,beta[n]);
    Dpol.resize(nphases, Dpol[n]);
    Spol.resize(nphases, Spol[n]);
    Jpol.resize(nphases, Jpol[n]);
    Dnem.resize(nphases, Dnem[n]);
    Snem.resize(nphases, Snem[n]);
    Jnem.resize(nphases, Jnem[n]);
    Knem.resize(nphases, Knem[n]);
    Wnem.resize(nphases, Wnem[n]);
    zetaQ.resize(nphases, zetaQ[n]);
    zetaS.resize(nphases, zetaS[n]);
    division.resize(nphases, division[n]);
    division_time.resize(nphases, division_time[n]);
    division_rate.resize(nphases, division_rate[n]);
    types.resize(nphases, types[n]);

    phi.resize(nphases, vector<double>(patch_N, 0.));
    phi_dx.resize(nphases, vector<double>(patch_N, 0.));
    phi_dy.resize(nphases, vector<double>(patch_N, 0.));
    phi_old.resize(nphases, vector<double>(patch_N, 0.));
    V.resize(nphases, vector<double>(patch_N, 0.));
    dphi.resize(nphases, vector<double>(patch_N, 0.));
    dphi_old.resize(nphases, vector<double>(patch_N, 0.));

    target_R.resize(nphases, target_R[n]);
    area.resize(nphases, 0.);
    overlap.resize(nphases, 0.);
    patch_min.resize(nphases, {0, 0});
    patch_max.resize(nphases, Size);
    com.resize(nphases, {0., 0.});
    com_prev.resize(nphases, {0., 0.});
    polarization.resize(nphases, {0., 0.});
    polarization_old.resize(nphases,{0.,0.});
    com_velocity.resize(nphases, {0., 0.});
    Fpassive.resize(nphases, {0., 0.});
    Fint.resize(nphases,{0.,0.});
    Frep.resize(nphases,{0.,0.});
    Fshape.resize(nphases, {0., 0.});
    Fnem.resize(nphases, {0., 0.});
    Fpol.resize(nphases, {0., 0.});
    com_x.resize(nphases, 0.);
    com_y.resize(nphases, 0.);
    S00.resize(nphases, 0.);
    S01.resize(nphases, 0.);
    Q00.resize(nphases, 0.);
    Q01.resize(nphases, 0.);
    offset.resize(nphases, {0u, 0u});
    theta_pol.resize(nphases, 0.);
    theta_pol_old.resize(nphases, 0.);
    theta_nem.resize(nphases, 0.);
    theta_nem_old.resize(nphases, 0.);
    vorticity.resize(nphases, 0.);
    tau.resize(nphases, 0.);
    division_counter.resize(nphases, 0.);
    sign_zetaQ.resize(nphases, 0.);
    sign_zetaS.resize(nphases, 0.);

    for(unsigned n=nphases-2; n<nphases; ++n)
    {
      a0.push_back(Pi*R[n]*R[n]);
      if(zetaQ[n]!=0.) sign_zetaQ[n] = zetaQ[n]>0. ? 1 : -1;
      if(zetaS[n]!=0.) sign_zetaS[n] = zetaS[n]>0. ? 1 : -1;
    }

    // create the cells
    const double theta0 = atan2(S01[n], S00[n])/2.;
    const double length = .5*sqrt(area[n]/Pi);
    const auto d = length * vec<double, 2> {cos(theta0), sin(theta0)};
    AddCell(nphases-2, coord(com[n] + d));
    AddCell(nphases-1, coord(com[n] - d));

    // -------------------------------------------------------------------------
    // Relax daughter cells

    vector<double> save_alpha(nphases, 0); swap(alpha, save_alpha);
    vector<double> save_beta(nphases, 0); swap(beta, save_beta);
    vector<double> save_zetaS(nphases, 0); swap(zetaS,  save_zetaS);
    vector<double> save_zetaQ(nphases, 0); swap(zetaQ,  save_zetaQ);
    vector<double> save_Dnem(nphases, 0); swap(Dnem,  save_Dnem);
    vector<double> save_Dpol(nphases, 0); swap(Dpol,  save_Dpol);
    vector<double> save_Jnem(nphases, 0); swap(Jnem,  save_Jnem);
    vector<double> save_Jpol(nphases, 0); swap(Jpol,  save_Jpol);
    vector<double> save_Knem(nphases, 0); swap(Knem,  save_Knem);
    vector<double> save_Wnem(nphases, 0); swap(Wnem,  save_Wnem);

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
    swap(beta, save_beta);
    swap(zetaQ, save_zetaQ);
    swap(zetaS, save_zetaS);
    swap(Jnem, save_Jnem);
    swap(Jpol, save_Jpol);
    swap(Dnem, save_Dnem);
    swap(Dpol, save_Dpol);
    swap(Knem, save_Knem);
    swap(Wnem, save_Wnem);

    cout << "END - phase " << n << "\n";

    // -------------------------------------------------------------------------
    // destroy the old fat mama

    // swap last and n
    SwapCells(nphases-1, n);

    // shrink memory
    nphases -= 1;

    gam.resize(nphases, gam[n]);
    mu.resize(nphases, mu[n]);
    R.resize(nphases, R[n]);
    xi.resize(nphases, xi[n]);
    alpha.resize(nphases, alpha[n]);
    beta.resize(nphases,beta[n]);
    Dpol.resize(nphases, Dpol[n]);
    Spol.resize(nphases, Spol[n]);
    Jpol.resize(nphases, Jpol[n]);
    Dnem.resize(nphases, Dnem[n]);
    Snem.resize(nphases, Snem[n]);
    Jnem.resize(nphases, Jnem[n]);
    Knem.resize(nphases, Knem[n]);
    Wnem.resize(nphases, Wnem[n]);
    zetaQ.resize(nphases, zetaQ[n]);
    zetaS.resize(nphases, zetaS[n]);
    division.resize(nphases, division[n]);
    division_time.resize(nphases, division_time[n]);
    division_rate.resize(nphases, division_rate[n]);
    types.resize(nphases, types[n]);

    phi.resize(nphases, vector<double>(patch_N, 0.));
    phi_dx.resize(nphases, vector<double>(patch_N, 0.));
    phi_dy.resize(nphases, vector<double>(patch_N, 0.));
    phi_old.resize(nphases, vector<double>(patch_N, 0.));
    V.resize(nphases, vector<double>(patch_N, 0.));
    dphi.resize(nphases, vector<double>(patch_N, 0.));
    dphi_old.resize(nphases, vector<double>(patch_N, 0.));

    target_R.resize(nphases);
    area.resize(nphases, 0.);
    patch_min.resize(nphases, {0, 0});
    patch_max.resize(nphases, Size);
    com.resize(nphases, {0., 0.});
    com_prev.resize(nphases, {0., 0.});
    polarization.resize(nphases, {0., 0.});
    polarization_old.resize(nphases,{0.,0.});
    com_velocity.resize(nphases, {0., 0.});
    Fpassive.resize(nphases, {0., 0.});
    Fint.resize(nphases,{0.,0.});
    Frep.resize(nphases,{0.,0.});
    Fshape.resize(nphases, {0., 0.});
    Fnem.resize(nphases, {0., 0.});
    Fpol.resize(nphases, {0., 0.});
    com_x.resize(nphases, 0.);
    com_y.resize(nphases, 0.);
    S00.resize(nphases, 0.);
    S01.resize(nphases, 0.);
    Q00.resize(nphases, 0.);
    Q01.resize(nphases, 0.);
    offset.resize(nphases, {0u, 0u});
    theta_pol.resize(nphases, 0.);
    theta_pol_old.resize(nphases, 0.);
    theta_nem.resize(nphases, 0.);
    theta_nem_old.resize(nphases, 0.);
    vorticity.resize(nphases, 0.);
    tau.resize(nphases, 0.);
    division_counter.resize(nphases, 0.);
    sign_zetaQ.resize(nphases, 0.);
    sign_zetaS.resize(nphases, 0.);
    a0.resize(nphases);
  }

  if(division_counter[n]<0)
  {
    R[n] = target_R[n]*(-division_counter[n])/double(nsubsteps*division_refract_time);

    if(unsigned(-division_counter[n]) == nsubsteps*division_refract_time)
      ResetDivisionCounter(n);
  }
}

void Model::ResetDivisionCounter(unsigned n)
{
  division_counter[n] = nsubsteps*(division_time[n] + 1
      + random_exponential(division_rate[n]));
}
