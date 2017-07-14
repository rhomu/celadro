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


