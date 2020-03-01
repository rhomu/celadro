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
#include "derivatives.hpp"

void Model::UpdateSubstrateDerivativesAtNode(unsigned k)
{
  const auto& s  = neighbors[k];  // stencil
  substrate_dxx_phi[k] = derivXX(substrate_phi, s);
  substrate_dyy_phi[k] = derivYY(substrate_phi, s);
}

void Model::UpdateSubstrateAtNode(unsigned k, bool store)
{
  const auto& s  = neighbors[k];  // stencil

  // fields and derivatives
  const auto& phi = substrate_phi[k];
  const auto& dxx_phi = substrate_dxx_phi[k];
  const auto& dyy_phi = substrate_dyy_phi[k];
  const auto& sxx = stress_xx[k];
  const auto& syy = stress_yy[k];
  const auto delta2_phi = laplacian(substrate_dxx_phi, s)
                          + laplacian(substrate_dyy_phi, s);

  // Swift-Hohenberg
  const auto dphi =
    + sxx*dxx_phi + syy*dyy_phi + substrate_gamma*delta2_phi
    - phi*(substrate_A + substrate_C*phi*phi);

  // predictor step only
  if(store)
  {
    substrate_phi_old[k] = phi;
    substrate_dphi_old[k] = dphi;
  }

  // predictor and corrector steps
  substrate_phi[k] = substrate_phi_old[k]
    + time_step*0.5*(dphi + substrate_dphi_old[k]);
}
