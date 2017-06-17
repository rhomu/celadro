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
#include "random.hpp"

using namespace std;

// truly random device to generate seed
random_device rd;
// pseudo random generator
// shlow
mt19937 gen(rd());
// super fasht
//ranlux24 gen(rd());
// uniform distribution on (0, 1)
uniform_real_distribution<> dist01(0, 1);
// uniform distribution on (0, 1)
uniform_real_distribution<> dist11(-1, 1);

// return random real, uniform distribution
double random_real(double min, double max)
{
  return uniform_real_distribution<>(min, max)(gen);
}

// return random real, gaussian distributed
double random_normal(double sigma)
{
  return normal_distribution<>(0., sigma)(gen);
}

// return geometric dist numbers, prob is p
unsigned random_geometric(double p)
{
  return geometric_distribution<>(p)(gen);
}

unsigned randu()
{
  return gen();
}

void set_seed(unsigned s)
{
  gen.seed(s);
}
