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

#ifndef WRITE_HPP_
#define WRITE_HPP_

#include "serialization.hpp"

/** Saves the complete frame of the system */
void WriteFrame(unsigned t);
/** Saves the complete frame of the system (old style)*/
void WriteFrame_old(unsigned t);
/** Write the simulation parameters */
void WriteParams();
/** Ask and maybe delete output files
 *
 * We do not overwrite any file! This is a mean of protection for light headed
 * grad students.
 * */
void ClearOutput();
/** Create temporary directory (in /tmp/ by default)
 *
 * The tmp directory used for output is /tmp/bunch-of-numbers, where bunch-
 * of-numbers is a bunch of numbers that has beens generated using a hash
 * function from the run name and a random salt.
 * */
void CreateOutputDir();

#endif//WRITE_HPP_
