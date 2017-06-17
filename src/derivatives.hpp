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

#ifndef DERIVATIVES_HPP_
#define DERIVATIVES_HPP_

#include <vector>
#include <array>

// =============================================================================
// Derivatives

/** Five-point finite difference derivative along the x direction */
inline double derivX(const std::vector<double>& arr,
                     const std::array<double, 9>& d,
                     const double stencil)
{
  return .5*(1-4*stencil)*(arr[d[1]]-arr[d[2]])
         + stencil*(arr[d[5]]-arr[d[6]]-arr[d[7]]+arr[d[8]]);
}

/** Five-point finite difference derivative along the x direction */
inline double derivY(const std::vector<double>& arr,
                     const std::array<double, 9>& d,
                     const double stencil)
{
  return .5*(1-4*stencil)*(arr[d[3]]-arr[d[4]])
         + stencil*(arr[d[5]]-arr[d[6]]+arr[d[7]]-arr[d[8]]);
}

/** Five-point finite difference second derivative along the x direction */
inline double derivXX(const std::vector<double>& arr,
                      const std::array<double, 9>& d,
                      const double stencil)
{
  return (1-4*stencil)*(arr[d[1]]+arr[d[2]])
   - 4*stencil*(arr[d[3]]+arr[d[4]])
   + 2*stencil*(arr[d[5]]+arr[d[6]]+arr[d[7]]+arr[d[8]])
   - 2*(1-4*stencil)*arr[d[0]];
}

/** Five-point finite difference second derivative along the z direction */
inline double derivYY(const std::vector<double>& arr,
                      const std::array<double, 9>& d,
                      const double stencil)
{
  return (1-4*stencil)*(arr[d[3]]+arr[d[4]])
   - 4*stencil*(arr[d[1]]+arr[d[2]])
   + 2*stencil*(arr[d[5]]+arr[d[6]]+arr[d[7]]+arr[d[8]])
   - 2*(1-4*stencil)*arr[d[0]];
}

/** Five-point finite difference derivative along x and z direction */
inline double derivXY(const std::vector<double>& arr,
                      const std::array<double, 9>& d,
                      const double stencil)
{
  return 0.5*(1-8*stencil)*(arr[d[1]]+arr[d[2]]+arr[d[3]]+arr[d[4]])
   + 2*stencil*(arr[d[5]]+arr[d[6]])
   - 0.5*(1-4*stencil)*(arr[d[7]]+arr[d[8]])
   - (1-8*stencil)*arr[d[0]];
}


/** Five-point finite difference laplacian (const) */
inline double laplacian(const std::vector<double>& arr,
                        const std::array<double, 9>& d,
                        const double stencil)
{
  return (1-4*stencil)*(arr[d[1]]+arr[d[2]]+arr[d[3]]+arr[d[4]])
         + 2*stencil*(arr[d[5]]+arr[d[6]]+arr[d[7]]+arr[d[8]])
         - 4*(1-2*stencil)*arr[d[0]];
}

/** Five-point finite difference flux
 *
 * This is just dx(ux*field) + dy(uy*field) with the same definition for the derivative as
 * derivX and derivY.
 * */
inline double flux(const std::vector<double>& arr,
                   const std::vector<double>& uxs,
                   const std::vector<double>& uzs,
                   const std::array<double, 9>& d,
                   const double stencil)
{
  return + .5*(1-4*stencil)*(uxs[d[1]]*arr[d[1]]-uxs[d[2]]*arr[d[2]])
         + stencil*(uxs[d[5]]*arr[d[5]]-uxs[d[6]]*arr[d[6]]-uxs[d[7]]*arr[d[7]]+uxs[d[8]]*arr[d[8]])
         + .5*(1-4*stencil)*(uzs[d[3]]*arr[d[3]]-uzs[d[4]]*arr[d[4]])
         + stencil*(uzs[d[5]]*arr[d[5]]-uzs[d[6]]*arr[d[6]]+uzs[d[7]]*arr[d[7]]-uzs[d[8]]*arr[d[8]]);
}

#endif//DERIVATIVES_HPP_
