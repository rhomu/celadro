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
#include "serialization.hpp"

using namespace std;

/** The number of spaces in indendation */
static const unsigned padding = 2;

oarchive::oarchive(std::ostream& stream_, string id, unsigned version)
  : stream(stream_)
{
  if(!stream.good()) throw bad_stream();
  // write initial brace
  open_group();
  // add description and version
  add("id", id);
  add("version", version);
  // open the data group
  add_key("data");
  open_group();
}

oarchive::~oarchive()
{
  // close the 'data' group
  close_group();
  // write final brace
  if(stream.good()) stream << endl << '}';
}

void oarchive::indent(unsigned n)
{
  level += n;
}

void oarchive::unindent(unsigned n)
{
  level -= n;
}

void oarchive::open_group(const std::string& obrace)
{
  // open brace
  stream << obrace;
  // indent one level
  indent();
  // next element is first
  first = true;
}

void oarchive::close_group(const string& cbrace)
{
  // unindent
  unindent();
  // new line
  stream << endl << string(padding*level, ' ');
  // close brace
  stream << cbrace;
  // not the first in the list
  first = false;
}

void oarchive::add_key(const string& key)
{
  new_line();
  stream << "\"" << key << "\" : ";
}

void oarchive::new_line()
{
  // add comma
  if(!first) stream << ',';
  else first = false;
  // new line
  stream << endl;
  // indent
  stream << string(padding*level, ' ');
}

template<>
void oarchive::add_element<const char*>(const char* const& t)
{
  std::stringstream ss;
  ss << "\"" << t << "\"";
  stream << ss.str();
}

template<>
void oarchive::add_element<std::string>(const std::string& t)
{
  add_element(t.c_str());
}
