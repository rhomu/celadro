/*
 * This file is part of CELADRO, Copyright (C) 2016-20, Romain Mueller
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
  sqrt_time_step = sqrt(time_step);

  // rectifies margin in case it is bigger than domain
  // and compensate for the boundary layer
  patch_margin = {
    min(margin, Size[0]/2 - 1 + (Size[0]%2)),
    min(margin, Size[1]/2 - 1 + (Size[1]%2))
  };
  // total size including bdry layer
  patch_size = 2u*patch_margin + 1u;
  patch_N = patch_size[0]*patch_size[1];

  // initialize memory for global fields
  walls.resize(N, 0.);
  walls_dx.resize(N, 0.);
  walls_dy.resize(N, 0.);
  walls_laplace.resize(N, 0.);
  sum.resize(N, 0.);
  pressure.resize(N, 0.);
  stress_xx.resize(N, 0.);
  stress_xy.resize(N, 0.);
  stress_yy.resize(N, 0.);
  sumA.resize(N, 0.);
  sumS00.resize(N, 0.);
  sumS01.resize(N, 0.);
  sumQ00.resize(N, 0.);
  sumQ01.resize(N, 0.);
  square.resize(N, 0.);
  thirdp.resize(N, 0.);
  fourthp.resize(N, 0.);
  P0.resize(N, 0.);
  P1.resize(N, 0.);
  U0.resize(N, 0.);
  U1.resize(N, 0.);

  // initialize substrate
  if(substrate_type)
  {
    substrate_phi.resize(N);
    substrate_phi_old.resize(N);
    substrate_dphi_old.resize(N);
    substrate_dxx_phi.resize(N);
    substrate_dyy_phi.resize(N);

    // set phi to small non-zero random values
    for(unsigned k=0; k<N; ++k)
      substrate_phi[k] = random_real(-1e-2, 1e-2);
  }

  // allocate memory for individual cells
  SetCellNumber(nphases);

  // ---------------------------------------------------------------------------

  // compute the sign of the activity coefficients
  if(zetaQ!=0.) sign_zetaQ = zetaQ>0. ? 1 : -1;
  if(zetaS!=0.) sign_zetaS = zetaS>0. ? 1 : -1;

  // compute tables
  for(unsigned i=0; i<Size[0]; ++i)
    com_x_table.push_back({ cos(-Pi+2.*Pi*i/Size[0]),
                            sin(-Pi+2.*Pi*i/Size[0]) });
  for(unsigned i=0; i<Size[1]; ++i)
    com_y_table.push_back({ cos(-Pi+2.*Pi*i/Size[1]),
                            sin(-Pi+2.*Pi*i/Size[1]) });

  // ---------------------------------------------------------------------------

  // check parameters
  for(unsigned n=0; n<nphases; ++n)
    if(margin<R) throw error_msg("Margin is too small, make it bigger than R.");

  // check birth boundaries
  if(birth_bdries.size()==0)
    birth_bdries = {0, Size[0], 0, Size[1]};
  else if(birth_bdries.size()!=4)
    throw error_msg("Birth boundaries have wrong format, see help.");

  if(wall_omega!=0)
    throw error_msg("Wall adhesion is not working for the moment.");

  if(substrate_type and BC!=0)
    throw error_msg("Substrate is only supported with PBC for the moment.");

  // ---------------------------------------------------------------------------

  a0 = Pi*R*R;
}

void Model::SetCellNumber(unsigned new_nphases)
{
  nphases = new_nphases;

  // allocate memory for qties defined on the patches
  phi.resize(nphases, vector<double>(patch_N, 0.));
  phi_dx.resize(nphases, vector<double>(patch_N, 0.));
  phi_dy.resize(nphases, vector<double>(patch_N, 0.));
  phi_old.resize(nphases, vector<double>(patch_N, 0.));
  V.resize(nphases, vector<double>(patch_N, 0.));
  dphi.resize(nphases, vector<double>(patch_N, 0.));
  dphi_old.resize(nphases, vector<double>(patch_N, 0.));

  // allocate memory for cell properties
  area.resize(nphases, 0.);
  patch_min.resize(nphases, {0, 0});
  patch_max.resize(nphases, Size);
  com.resize(nphases, {0., 0.});
  com_prev.resize(nphases, {0., 0.});
  polarization.resize(nphases, {0., 0.});
  velocity.resize(nphases, {0., 0.});
  Fpressure.resize(nphases, {0., 0.});
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
  delta_theta_pol.resize(nphases, 0.);
  theta_nem.resize(nphases, 0.);
  theta_nem_old.resize(nphases, 0.);
  vorticity.resize(nphases, 0.);
  tau.resize(nphases, 0.);
}

void Model::InitializeNeighbors()
{
  neighbors.resize(N);
  neighbors_patch.resize(patch_N);

  // define the neighbours, accounting for the periodic boundaries
  for(unsigned k=0; k<N; ++k)
  {
    const unsigned x = GetXPosition(k);
    const unsigned y = GetYPosition(k);
    for(int dx=-1; dx<=1; ++dx)
      for(int dy=-1; dy<=1; ++dy)
        neighbors[k][dx][dy] = GetIndex({ (x+Size[0]+dx)%Size[0],
                                          (y+Size[1]+dy)%Size[1] });
  }

  // define the neighbours, accounting for the boundary layer
  const unsigned lx = patch_size[0];
  const unsigned ly = patch_size[1];
  for(unsigned k=0; k<patch_N; ++k)
  {
    const unsigned x = k/ly;
    const unsigned y = k%ly;
    for(int dx=-1; dx<=1; ++dx)
    {
      for(int dy=-1; dy<=1; ++dy)
      {
        const unsigned u = (x+lx+dx)%lx;
        const unsigned v = (y+ly+dy)%ly;
        neighbors_patch[k][dx][dy] = v + u*ly;
      }
    }
  }
}
