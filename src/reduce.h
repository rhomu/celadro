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

#ifndef REDUCE_H_
#define REDUCE_H_

#include "cuda.h"

// =============================================================================
// Warp reduce, from:
// https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/

template<class T>
__inline__ __device__
void warpReduceSum(T& val)
{
  for (int offset = WarpSize/2; offset > 0; offset /= 2)
    val += __shfl_down(val, offset);
}

template<class T, size_t D>
__inline__ __device__
void warpReduceSum(vec<T, D> &val)
{
  for (int offset = WarpSize/2; offset > 0; offset /= 2)
    for(int i=0; i<D; ++i)
      val[i] += __shfl_down(val[i], offset);
}

template<class T>
__inline__ __device__
void blockReduceSum(T& val)
{
  static __shared__ T shared[WarpSize];
  int lane = threadIdx.x % WarpSize;
  int wid  = threadIdx.x / WarpSize;

  warpReduceSum(val);     // Each warp performs partial reduction

  if (lane==0) shared[wid]=val; // Write reduced value to shared memory

  __syncthreads();              // Wait for all partial reductions

  // read from shared memory only if that warp existed
  if(threadIdx.x < blockDim.x / warpSize)
    val = shared[lane];
  else 
    val = 0;

  if (wid==0) warpReduceSum(val); // Final reduce within first warp
}

#endif//REDUCE_H_
