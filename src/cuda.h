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

/** There are two different compiler flags associated with CUDA:
 *
 * _CUDA is defined when compiling device source code (using NVCC)
 * _CUDA_ENABLED is defined project-wise when cuda is enabled, i.e. also for
 *               host source code.
 **/

#ifndef CUDA_HPP_
#define CUDA_HPP_

#ifdef _CUDA

#define _CUDA_ENABLED

#ifdef _OPENMP
#error "Cuda can not be used along with OpenMP"
#endif

#define WarpSize 32

#define ThreadsPerBlock 1024

#define CUDA_host_device __host__ __device__

#else

#define CUDA_host_device 

#endif
#endif//CUDA_HPP_
