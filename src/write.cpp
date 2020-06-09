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
#include "model.hpp"
#include "files.hpp"

using namespace std;

void Model::WriteFrame(unsigned t)
{
  // construct output name
  const string oname = inline_str(output_dir, "frame", t, ".json");

  // write
  {
    stringstream buffer;
    {
      oarchive ar(buffer, "frame", 1);
      // serialize
      SerializeFrame(ar);

      if(ar.bad_value()) throw error_msg("nan found while writing file.");
    }
    // dump to file
    std::ofstream ofs(oname.c_str(), ios::out);
    ofs << buffer.rdbuf();
  }

  // compress
  if(compress) compress_file(oname, oname);
}

void Model::WriteParams()
{
  // a name that makes sense
  const string oname = inline_str(output_dir, "parameters.json");

  // write
  {
    stringstream buffer;
    {
      // serialize
      oarchive ar(buffer, "parameters", 1);
      // ...program parameters...
      ar & auto_name(Size)
         & auto_name(BC)
         & auto_name(nsteps)
         & auto_name(nsubsteps)
         & auto_name(ninfo)
         & auto_name(nstart);
      // ...and model parameters
      SerializeParameters(ar);

      if(ar.bad_value()) throw error_msg("nan found while writing file.");
    }
    // dump to file
    std::ofstream ofs(oname.c_str(), ios::out);
    ofs << buffer.rdbuf();
  }

  // compress
  if(compress) compress_file(oname, oname);
}

void Model::ClearOutput()
{
  output_dir = inputname + ( inputname.back()=='/' ? "data/" : "/data/" );

  // extension of single files
  string ext = compress ? ".json.zip" : ".json";

  if(fs::is_directory(output_dir))
  {
    if(not fs::is_empty(output_dir) and not force_delete)
    {
      auto ans = ask(" remove output files in directory '" + output_dir + "'?");

      if(not ans)
        throw error_msg("output files already exist in directory '",
            output_dir, "'.");
    }

    fs::remove_all(output_dir);
  }

  fs::create_directory(output_dir);
}
void Model::WriteConfig(){
  // a name that makes sense
  const string oname = inline_str(output_dir, "init_config0.json");

  // write
  {
    stringstream buffer;
    {
      // serialize
      oarchive ar(buffer, "init_config", 1);
      // ...program parameters...
      ar & auto_name(Size)
         & auto_name(BC);
      SerializeParameters(ar);

      if(ar.bad_value()) throw error_msg("nan found while writing file.");
    }
    // dump to file
    std::ofstream ofs(oname.c_str(), ios::out);
    ofs << buffer.rdbuf();
  }
}

