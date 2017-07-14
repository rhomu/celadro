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

void Model::Initialize()
{
  N = Size[0]*Size[1];

  // initialize memory
  walls.resize(N, 0.);
  walls_dx.resize(N, 0.);
  walls_dy.resize(N, 0.);
  walls_laplace.resize(N, 0.);
  sum.resize(N, 0.);
  sum_cnt.resize(N, 0.);
  square.resize(N, 0.);
  square_cnt.resize(N, 0.);
  Px.resize(N, 0.);
  Py.resize(N, 0.);
  Theta.resize(N, 0.);
  Q00.resize(N, 0.);
  Q01.resize(N, 0.);
  Px_cnt.resize(N, 0.);
  Py_cnt.resize(N, 0.);
  Theta_cnt.resize(N, 0.);
  Q00_cnt.resize(N, 0.);
  Q01_cnt.resize(N, 0.);
  P.resize(N, 0.);
  P_cnt.resize(N, 0.);

  // the patch is the region over which we compute each cell
  const unsigned PatchSize = N;//2*margin+1;

  phi.resize(nphases, vector<double>(PatchSize, 0.));
  phi_old.resize(nphases, vector<double>(PatchSize, 0.));
  V.resize(nphases, vector<double>(PatchSize, 0.));
  potential.resize(nphases, vector<double>(PatchSize, 0.));
  potential_old.resize(nphases, vector<double>(PatchSize, 0.));

  area.resize(nphases, 0.);
  area_cnt.resize(nphases, 0.);
  domain_min.resize(nphases, {0, 0});
  domain_max.resize(nphases, Size);
  com.resize(nphases, {0., 0.});
  com_prev.resize(nphases, {0., 0.});
  pol.resize(nphases, {0., 0.});
  velp.resize(nphases, {0., 0.});
  velc.resize(nphases, {0., 0.});
  velf.resize(nphases, {0., 0.});
  vel.resize(nphases, {0., 0.});
  com_x.resize(nphases, 0.);
  com_y.resize(nphases, 0.);
  c.resize(nphases, 0);
  S00.resize(nphases, 0.);
  S01.resize(nphases, 0.);
  S_order.resize(nphases, 0.);
  S_angle.resize(nphases, 0.);
  theta.resize(nphases, 0.);

  // pre-compute coefficients
  C1 = 60./lambda/lambda;
  C2 = Pi*R*R;
  C3 = C1/lambda/lambda;

  // extend the parameters with the last given value
  gam.resize(nphases, gam.back());
  mu.resize(nphases, mu.back());

  // compute tables
  for(unsigned i=0; i<Size[0]; ++i)
    com_x_table.push_back({ cos(-Pi+2.*Pi*i/Size[0]), sin(-Pi+2.*Pi*i/Size[0]) });
  for(unsigned i=0; i<Size[1]; ++i)
    com_y_table.push_back({ cos(-Pi+2.*Pi*i/Size[1]), sin(-Pi+2.*Pi*i/Size[1]) });

  // check birth boundaries
  if(birth_bdries.size()==0)
    birth_bdries = {0, Size[0], 0, Size[1]};
  else if(birth_bdries.size()!=4)
    throw error_msg("Birth boundaries have wrong format, see help.");

  // check margin size
  if(2*margin+1>Size[0] or 2*margin+1>Size[1])
    throw error_msg("Margin is too large for domain size.");

  // set flags
  division = (division_rate!=0.);
}

void Model::InitializeNeighbors()
{
  neighbors.resize(N);

  // define the neighbours, accounting for the periodic boundaries
  for(unsigned k=0; k<N; ++k)
  {
    const unsigned x = GetXPosition(k);
    const unsigned y = GetYPosition(k);
    for(int dx=-1; dx<=1; ++dx)
      for(int dy=-1; dy<=1; ++dy)
        neighbors[k][dx][dy] = GetDomainIndex( (x+Size[0]+dx)%Size[0], (y+Size[1]+dy)%Size[1] );
  }
}
