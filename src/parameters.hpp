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

#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

/** Declare and parse all program options
  *
  * This function takes care of the option parsing. This tend to be a little
  * intricate due to the different ways the user can interact with both the
  * command line and the runcard file. Moreover, we have the complication that
  * the different models generally defines options that can only be parsed after
  * the current model is known.
  * */
void ParseProgramOptions(int ac, char **av);

/** Process program options
 *
 * This functions makes small adjustments and checks to the inputed program
 * options obtained from ParseProgramOptions and is called just afterwards. Note
 * that most options do not need any processing.
 * */
void ProcessProgramOptions();

/** Print simulation parameters
  *
  * Simply prints all set parameters with a (somewhat) nice formatting. */
void PrintProgramOptions();

#endif//OPTIONS_HPP_
