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

#ifndef RANDOM_HPP_
#define RANDOM_HPP_

/** Return random real between min and max, uniform distribution */
double random_real(double min=0., double max=1.);

/** Return random real, normal distribution with variance sigma and zero mean */
double random_normal(double sigma=1.);

/** Return geometric distributed integers */
unsigned random_geometric(double p);

/** Return random iunsigned int */
unsigned randu();

/** Set seed of the random number generator */
void set_seed(unsigned);

#endif//RANDOM_HPP_
