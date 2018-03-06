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

#ifndef VEC_HPP_
#define VEC_HPP_

#include <cmath>
#include <array>
#include <iostream>

/** Simple euclidean vector
 *
 * The all-time most famous class among C++ students! This is simply an array
 * with euclidean operations built on top.
 * */
template<class T, size_t D>
struct vec
{
  /** Individual components */
  std::array<T, D> data;

  vec() = default;
  vec(const vec& v) = default;
  vec& operator=(const vec& v) = default;

  vec& operator+=(const vec& v)
  {
    for(size_t i = 0; i<D; ++i)
      data[i] += v.data[i];
    return *this;
  }

  vec& operator-=(const vec& v)
  {
    for(size_t i = 0; i<D; ++i)
      data[i] -= v.data[i];
    return *this;
  }

  vec& operator*=(const T& t)
  {
    for(auto& e : data)
      e *= t;
    return *this;
  }

  vec& operator/=(const T& t)
  {
    for(auto& e : data)
      e /= t;
    return *this;
  }

  bool operator!=(const vec& v) const
  {
    for(size_t i = 0; i<D; ++i)
      if(data[i] != v.data[i])
        return true;
    return false;
  }

  bool operator==(const vec& v) const
  { return not (*this!=v); }

  vec operator+(const vec& v) const
  {
    vec ret;
    for(size_t i = 0; i<D; ++i)
      ret[i] = data[i] + v.data[i];
    return ret;
  }

  vec operator-(const vec& v) const
  {
    vec ret;
    for(size_t i = 0; i<D; ++i)
      ret[i] = data[i] - v.data[i];
    return ret;
  }

  T operator*(const vec& v) const
  {
    T ret {0};
    for(size_t i = 0; i<D; ++i)
      ret += data[i]*v.data[i];
    return ret;
  }

  T& operator[](size_t i)
  { return data[i]; }

  const T& operator[](size_t i) const
  { return data[i]; }

  /** Component-wise comparison operator */
  bool operator>(const vec& v) const
  {
    for(size_t i=0; i<D; ++i)
      if(data[i]<=v.data[i])
        return false;
    return true;
  }

  /** Component-wise comparison operator */
  bool operator<(const vec& v) const
  {
    for(size_t i=0; i<D; ++i)
      if(data[i]>=v.data[i])
        return false;
    return true;
  }

  /** Component-wise comparison operator */
  bool operator>=(const vec& v) const
  {
    for(size_t i=0; i<D; ++i)
      if(data[i]<v.data[i])
        return false;
    return true;
  }

  /** Component-wise comparison operator */
  bool operator<=(const vec& v) const
  {
    for(size_t i=0; i<D; ++i)
      if(data[i]>v.data[i])
        return false;
    return true;
  }

  /** Square */
  T sq() const
  { return (*this)*(*this); }

  /** Absolute value */
  T abs() const
  { return sqrt(abs(sq())); }

  /** Return unit vector */
  vec unit_vector() const
  { return *this/sqrt(sq()); }

  template<class U, size_t E> friend
  vec<U, E> operator+(const U&, const vec<U, E>&);
  template<class U, size_t E> friend
  vec<U, E> operator+(const vec<U, E>&, const U&);
  template<class U, size_t E> friend
  vec<U, E> operator-(const U&, const vec<U, E>&);
  template<class U, size_t E> friend
  vec<U, E> operator-(const vec<U, E>&, const U&);
  template<class U, size_t E> friend
  vec<U, E> operator*(const U&, const vec<U, E>&);
  template<class U, size_t E> friend
  vec<U, E> operator*(const vec<U, E>&, const U&);
  template<class U, size_t E> friend
  vec<U, E> operator/(const vec<U, E>&, const U&);
  template<class U, size_t E> friend
  std::ostream& operator<<(std::ostream&, const vec<U, E>&);

  template<class U>
  explicit operator vec<U, D>() const
  {
    vec<U, D> ret;
    for(size_t i = 0; i<D; ++i)
      ret[i] = U(data[i]);
    return ret;
  }

  // for serialization
  using value_type = T;
  typename std::array<T, D>::iterator begin() { return data.begin(); }
  typename std::array<T, D>::iterator end()   { return data.end(); }
  typename std::array<T, D>::const_iterator begin() const { return data.begin(); }
  typename std::array<T, D>::const_iterator end()   const { return data.end(); }
  size_t size() const { return D; }
};

// =============================================================================
// External functions implementation

template<class T, size_t D>
vec<T, D> operator+(const vec<T, D>& v, const T& t)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i] + t;
  return ret;
}

template<class T, size_t D>
vec<T, D> operator+(const T& t, const vec<T, D>& v)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i] + t;
  return ret;
}

template<class T, size_t D>
vec<T, D> operator-(const vec<T, D>& v, const T& t)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i] - t;
  return ret;
}

template<class T, size_t D>
vec<T, D> operator-(const T& t, const vec<T, D>& v)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i]-t;
  return ret;
}

template<class T, size_t D>
vec<T, D> operator*(const vec<T, D>& v, const T& t)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i]*t;
  return ret;
}

template<class T, size_t D>
vec<T, D> operator*(const T& t, const vec<T, D>& v)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i]*t;
  return ret;
}

template<class T, size_t D>
vec<T, D> operator/(const vec<T, D>& v, const T& t)
{
  vec<T, D> ret;
  for(size_t i = 0; i<D; ++i)
    ret[i] = v.data[i]/t;
  return ret;
}

/** Ostream output
  *
  * Prints all components in this form: [0, 1, ...]. */
template<class T, size_t D>
std::ostream& operator<<(std::ostream& stream, const vec<T, D>& v)
{
  stream << '[';

  if(D>0)
  {
    stream << v[0];
    for(size_t i = 1; i<D; ++i)
      stream << ", " << v[i];
  }

  return stream << ']';
}

/** Modulo operator for vec<unisgned, D>
 *
 * Applies modulo component-wise.
 * */
template<size_t D>
vec<unsigned, D> operator%(const vec<unsigned, D>& a, const vec<unsigned, D>& b)
{
  auto ret = a;
  for(size_t i = 0; i<D; ++i)
    ret.data[i] %= b.data[i];
  return ret;
}

/** Wrap vector around periodic boundaries */
template<class T, size_t D>
inline vec<T, D> wrap(vec<T, D> v, const vec<T, D>& L)
{
  for(size_t i=0; i<D; ++i)
    v[i] = wrap(v[i], L[i]);
  return v;
}

#endif//VEC_HPP_
