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

#include "stencil.hpp"

// =============================================================================
// Derivatives

/** Symmetric finite difference derivative along the x direction */
inline double derivX(const field& f, const stencil& s)
{
  return .5*( f[s[+1][0]] - f[s[-1][0]] );
}

/** Symmetric finite difference derivative along the y direction */
inline double derivY(const field& f, const stencil& s)
{
  return .5*( f[s[0][+1]] - f[s[0][-1]] );
}

/** Five-point finite difference laplacian */
inline double laplacian(const field& f, const stencil& s)
{
  return f[s[+1][0]] + f[s[0][+1]] + f[s[-1][0]] + f[s[0][-1]] - 4.*f[s[0][0]];
}

/** Symmetric finite difference derivative along the x direction (arrays) */
inline double derivX(double *f, const stencil& s)
{
  return .5*( f[s[+1][0]] - f[s[-1][0]] );
}

/** Symmetric finite difference derivative along the y direction (arrays) */
inline double derivY(double *f, const stencil& s)
{
  return .5*( f[s[0][+1]] - f[s[0][-1]] );
}

/** Five-point finite difference laplacian (arrays) */
inline double laplacian(double *f, const stencil& s)
{
  return f[s[+1][0]] + f[s[0][+1]] + f[s[-1][0]] + f[s[0][-1]] - 4.*f[s[0][0]];
}


#endif//DERIVATIVES_HPP_
