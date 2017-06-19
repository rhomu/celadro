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

#include "threads.hpp"

extern unsigned nthreads;

/** Init multi-threading
  *
  * Somewhow omp_get_num_threads() is not working properly with gcc... so we
  * need to use this trick to get the standard number of threads. Or I am too
  * dumb to use omp...
  * */
void SetThreads()
{
  // if nthreads is 1 we use the default number of threads from OpenMP
  if(nthreads == 1)
  {
    // count the number of OpenMP threads
    unsigned count = 0;
    PRAGMA_OMP(omp parallel)
    {
      PRAGMA_OMP(omp atomic)
      ++count;
    }
    nthreads = count;
  }
}
