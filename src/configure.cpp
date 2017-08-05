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

using namespace std;

void Model::AddCellAtNode(unsigned n, unsigned k, const coord& center)
{
  const unsigned xk = GetXPosition(k);
  const unsigned yk = GetYPosition(k);
  const unsigned q  = GetPhiIndex(n, xk, yk);

  // we create smaller cells that will then relax
  // this improves greatly the stability at the first steps
  const auto radius = max(R/2., 4.);

  // round shape (do not wrap if no PBC)
  if(
      (BC==0 and pow(wrap(diff(yk, center[1]), Size[1]), 2)
       + pow(wrap(diff(xk, center[0]), Size[0]), 2)<=ceil(radius*radius))
      or
      (BC>=1 and pow(diff(yk, center[1]), 2)
       + pow(diff(xk, center[0]), 2)<=ceil(radius*radius))
    )
  {
    phi[n][q]     = 1.;
    phi_old[n][q] = 1.;
    area[n]      += 1.;
    square[k]    += 1.;
    sum[k]       += 1.;
  }
  else
  {
    phi[n][q]     = 0.;
    phi_old[n][q] = 0.;
  }
}

void Model::AddCell(unsigned n, const coord& center)
{
  // update patch coordinates
  patch_min[n] = (center+Size-patch_size)%Size;
  patch_max[n] = (center+patch_size+1u)%Size;

  // create the cells at the centers we just computed
  UpdateDomain(&Model::AddCellAtNode, n ,center);

  c[n]     = c0;
  theta[n] = 2.*Pi*random_real();
  com[n]   = vec<double, 2>(center);
}

void Model::Configure()
{
  // ===========================================================================
  // adding cells at random while trying to keep their center non overlapping
  if(init_config=="random" and BC==0)
  {
    // target radius for spacing between cells
    unsigned radius = sqrt(double(N/nphases)/Pi);
    // list of all centers
    vector<coord> centers;

    for(unsigned n=0; n<nphases; ++n)
    {
      // generate new center while trying to keep safe distance
      while(true)
      {
        const coord center = {
          static_cast<unsigned>(birth_bdries[0]
              +random_real()*(birth_bdries[1]-birth_bdries[0])),
          static_cast<unsigned>(birth_bdries[2]
              +random_real()*(birth_bdries[3]-birth_bdries[2]))
        };

        bool is_overlapping = false;
        for(const auto& c : centers)
        {
          if(pow(wrap(diff(center[0], c[0]), Size[0]), 2)
              + pow(wrap(diff(center[1], c[1]), Size[1]), 2) < 0.9*radius*radius)
          {
            is_overlapping = true;
            break;
          }
        }

        if(!is_overlapping)
        {
          centers.emplace_back(center);
          break;
        }
      }

      // add cell
      AddCell(n, centers.back());
    }
  }
  // ===========================================================================
  // same but with walls: we need to be careful not to create cells on the wall
  else if(init_config=="random" and BC>=1)
  {
    // target radius for spacing between cells
    unsigned radius = sqrt(double(N/nphases)/Pi);
    // list of all centers
    vector<coord> centers;

    for(unsigned n=0; n<nphases; ++n)
    {
      // generate new center while trying to keep safe distance
      while(true)
      {
        const coord center = {
          static_cast<unsigned>(birth_bdries[0]
              +random_real()*(birth_bdries[1]-birth_bdries[0])),
          static_cast<unsigned>(birth_bdries[2]
              +random_real()*(birth_bdries[3]-birth_bdries[2]))
        };

        // detect walls
        // ... box only
        if(BC==1)
          if(center[0]<0.9*R or Size[0]-center[0]<0.9*R) continue;
        // ... box and channel
        if(BC==1 or BC==2)
          if(center[1]<0.9*R or Size[1]-center[1]<0.9*R) continue;
        // ... ellipse
        if(BC==3)
        {
          // compute distance from the elliptic wall
          // ... angle of the current point (from center of the patch)
          const auto theta = atan2(Size[1]/2.-center[1], Size[0]/2.-center[0]);
          // ... small helper function to compute radius
          auto rad = [](auto x, auto y) { return sqrt(x*x + y*y); };
          // ... distance is the difference between wall and current point
          const auto d = rad(Size[0]/2.*cos(theta), Size[1]/2.*sin(theta))
                        -rad(Size[0]/2.-center[0], Size[1]/2.-center[1]);

          if(d<0.9*R) continue;
        }

        // overlapp between cells
        bool is_overlapping = false;
        for(const auto& c : centers)
        {
          if(pow(wrap(diff(center[0], c[0]), Size[0]), 2)
              + pow(wrap(diff(center[1], c[1]), Size[1]), 2) < 0.9*radius*radius)
          {
            is_overlapping = true;
            break;
          }
        }

        if(!is_overlapping)
        {
          centers.emplace_back(center);
          break;
        }
      }

      // add cell
      AddCell(n, centers.back());
    }
  }
  // ===========================================================================
  // cluster of close cells in the center
  else if(init_config=="cluster")
  {
    const double theta  = 2*Pi/nphases;
    const double radius = R + nphases - 2;

    for(unsigned n=0; n<nphases; ++n)
      AddCell(n, {unsigned(Size[0]/2+radius*(cos(n*theta)+noise*random_real())),
                  unsigned(Size[1]/2+radius*(sin(n*theta)+noise*random_real())) });
  }
  // ===========================================================================
  // single cell in the middle
  else if(init_config=="single")
  {
    if(nphases!=1)
      throw error_msg("error: initial conditions require "
                      "nphases=1.");

    AddCell(0, {Size[0]/2, Size[1]/2});
  }
  else throw error_msg("error: initial configuration '",
      init_config, "' unknown.");
}

void Model::ConfigureWalls()
{
  switch(BC)
  {
  case 0:
    // no walls (pbc)
    for(unsigned k=0; k<N; ++k)
      walls[k] = 0;
    break;
  case 1:
    // Exponentially falling phase-field:
    for(unsigned k=0; k<N; ++k)
    {
      const double x = GetXPosition(k);
      const double y = GetYPosition(k);

      // this is the easiest way: each wall contributes as an exponentially
      // falling potential and we do not care about overalps
      walls[k] =   exp(-y/wall_thickness)
                 + exp(-x/wall_thickness)
                 + exp(-(Size[0]-1-x)/wall_thickness)
                 + exp(-(Size[1]-1-y)/wall_thickness);
    }
    break;
  // Same as above but channel.
  case 2:
    for(unsigned k=0; k<N; ++k)
    {
      const auto y = GetYPosition(k);

      // exponentially falling on both sides
      walls[k] = exp(-double(y)/wall_thickness)
        + exp(-double(Size[1]-y-1)/wall_thickness);
    }
    break;
  // ellipse!
  case 3:
    for(unsigned k=0; k<N; ++k)
    {
      const auto x = GetXPosition(k);
      const auto y = GetYPosition(k);

      // compute distance from the elliptic wall
      // ... angle of the current point (from center of the domain)
      const auto theta = atan2(Size[1]/2.-y, Size[0]/2.-x);
      // ... small helper function to compute radius
      const auto rad = [](auto x, auto y) { return sqrt(x*x + y*y); };
      // ... distance is the difference between wall and current point
      const auto d = rad(Size[0]/2.*cos(theta), Size[1]/2.*sin(theta))
                    -rad(Size[0]/2.-x, Size[1]/2.-y);
      // set the wall
      if(d<0)
        walls[k] = 1.;
      else
        walls[k] = exp(-d/wall_thickness);
    }
    break;
  default:
    throw error_msg("boundary condition unknown.");
  }

  // pre-compute derivatives
  for(unsigned k=0; k<N; ++k)
  {
    const auto& s = neighbors[k];

    walls_dx[k] = derivX(walls, s);
    walls_dy[k] = derivY(walls, s);
    walls_laplace[k] = laplacian(walls, s);
  }
}

