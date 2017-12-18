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
#include "files.hpp"
using namespace std;

void create_directory(const string& dir)
{
  const int ret = system(inline_str("mkdir -p ", dir).c_str());
  if(ret) throw error_msg("can not create output directory, mkdir returned ", ret, ".");
}

void remove_file(const string& fname)
{
  const int ret = system(inline_str("rm -rf ", fname).c_str());
  if(ret) throw error_msg("rm returned non-zero value ", ret, ".");
}

void compress_file(const string& iname, const string& oname)
{
  const int ret = system(inline_str("zip -jm ", oname, ".zip ", iname,
                                    " > /dev/null 2>&1").c_str());
  if(ret!=0) throw error_msg("zip non-zero return value ", ret, ".");
}
