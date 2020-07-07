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

#ifndef TOOLS_HPP_
#define TOOLS_HPP_

#include <cmath>
#include <sstream>
#include <iostream>
#include <limits>
#include <vector>

/** Modulo function that works correctly with negative values */
template<class T>
inline T modu(const T& num, const T& div)
{
  if(num < 0) return div + num%div;
  else return num%div;
}

/** Modulo function that works correctly with negative values */
inline double modu(double num, double div)
{
  if(num < 0) return div + std::fmod(num, div);
  else return std::fmod(num, div);
}

/** Unsigned difference
 *
 * This function can be used with unsigned types with no fear of loosing the
 * sign.
 * */
template<class T>
inline T diff(const T& a, const T& b)
{
  return a>b ? a-b : b-a;
}

/** Check that two doubles are equal
 *
 * This somewhat complicated expression avoids false negative from finite
 * precision of the floating point arithmetic. Precision threshold can be
 * set using the error_factor parameter.
 * */
inline bool check_equal(double a, double b, double error_factor=1.)
{
  return a==b ||
    std::abs(a-b)<std::abs(std::min(a,b))*std::numeric_limits<double>::epsilon()*
                  error_factor;
}

/** Wrap point around periodic boundaries
 *
 * This function is useful when dealing with periodic boundary conditions,
 * simply returns the minimum of x%L and L-x%L
 * */
template<class T, class U>
inline T wrap(const T& x, const U& L)
{
  return std::min(modu(x, L), L-modu(x, L));
}

/** Specialization for double */
inline double wrap(double x, double L)
{
  const auto y = modu(x, L);
  if(std::abs(y)<std::abs(L-y))
    return y;
  else
    return L-y;
}

/** Set if smaller
 *
 * This function sets the first variable to be equal to the second if it is
 * smaller.
 * */
template<class T, class U>
inline void set_if_smaller(T& dst, const U& src)
{
  if(src<dst) dst = src;
}

/** Set if bigger
 *
 * This function sets the first variable to be equal to the second if it is
 * bigger.
 * */
template<class T, class U>
inline void set_if_bigger(T& dst, const U& src)
{
  if(src>dst) dst = src;
}

namespace detail
{
  /** Convert to strig and catenate arguments */
  template<class Head>
  void inline_str_add_args(std::ostream& stream, Head&& head)
  {
    stream << std::forward<Head>(head);
  }
  /** Convert to strig and catenate arguments */
  template<class Head, class... Tail>
  void inline_str_add_args(std::ostream& stream, Head&& head, Tail&&... tail)
  {
    stream << std::forward<Head>(head);
    inline_str_add_args(stream, std::forward<Tail>(tail)...);
  }
} // namespace detail

/** Convert any number of arguments to string and catenate
 *
 * It does pretty much what is advertised. Look at the code if you want to learn
 * some pretty neat modern C++.
 * */
template<class... Args>
std::string inline_str(Args&&... args)
{
  std::stringstream s;
  detail::inline_str_add_args(s, std::forward<Args>(args)...);
  return s.str();
}

/** Convert iterable to string of the form {a,b,c,...} */
template<class T>
std::string vec2str(const T& iterable)
{
  std::stringstream s;
  s << '{';
  for(auto it = begin(iterable);;)
  {
    s << *it;
    if(++it==end(iterable)) break;
    s << ',';
  }
  s << '}';

  return s.str();
}

/** Split string
 *
 * My most famous contribution to stack overflow.
 * */
inline std::vector<std::string> split(const std::string& s, char sep=' ')
{
  std::vector<std::string> words;
  for(size_t p=0, q=0; p!=s.npos; p=q)
    words.push_back(s.substr(p+(p!=0), (q=s.find(sep, p+1))-p-(p!=0)));
  return words;
}

#endif//TOOLS_HPP_
