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

/** Small json serialization library.
  *
  * Minimal json library for arbitrary types. It goes together with the python
  * file archive.py that allows to import automatically the data dumped using
  * this library.
  *
  * Usage is similar to boost::serialization, but the API is minimal. The
  * implementation is also inspired from this library. The main choice is to
  * avoid RTI and use templates. This is more efficient, but comes with a lot
  * of difficulties because template functions can not be virtual in C++. Hence
  * we have to resort to dirty template tricks (and admit it you like dirty
  * tricks). This is a very good opportunity to write code that is: (i) overly
  * complicated and (ii) unmaintainable by any standard.
  *
  * Remark: We can not deduce names automatically because these are
  * implementation defined (even through <typeinfo>). This would mean that an
  * output file produced on a computer might not be usable on another, which
  * is uncool.
  *
  * TODO:
  * * Declare string types
  * * InputArchive
  *
  * Be careful, part of this code is serious C++.
  *
  * Romain Mueller, name dot surname at gmail dot com
  */

#ifndef SERIALIZATION_HPP_
#define SERIALIZATION_HPP_

// =============================================================================
// Tools

/** Bad stream exception
  *
  * In rare cases the stream might smell funny before the archive destructor is
  * called. We throw an exception instead of creating a seg fault.
  * */
struct bad_stream : public std::exception
{
  const char * what () const throw ()
  { return "can not use stream"; }
};

/** Pairs a variable with its display name
  *
  * We need to provide names in json! We do this using std::pair. The name of a
  * variable can be anything but if you want the same name as in the code you
  * can use this macro.
  * */
#define auto_name(obj) std::pair<decltype(obj)&, std::string> {obj, #obj}

// =============================================================================
// Type traits (using SFINAE)

namespace detail
{
  /** Check if type is iterable
    *
    * We use SFINAE to check at compile time if a type is iterable, i.e. if it
    * implements begin, end, !=, operator++ and dereferencement for iterators.
    * */
  template<typename T>
  auto is_iterable_impl(int)
    -> decltype (
        // begin/end and operator !=
        std::begin(std::declval<T&>()) != std::end(std::declval<T&>()),
        // operator ++ on begin iterator
        ++std::declval<decltype(std::begin(std::declval<T&>()))&>(),
        // operator* on begin iterator
        *std::begin(std::declval<T&>()),
        // value type
        std::declval<typename T::value_type>(),
        std::true_type {}
      );

  template<typename T>
  std::false_type is_iterable_impl(...);

  /** Check if type has a serialize() function
    *
    * We use SFINAE to detect at compile time if a certain type implements a
    * serialize function.
    * */
  template<class T>
  auto is_serializable_impl(int)
    -> decltype(
        std::declval<T>().serialize(std::declval<int&>()),
        std::true_type {}
      );

  template<class>
  std::false_type is_serializable_impl(...);

  /** Wrapper to serialization
    *
    * This objects allows us to differentiate objects that implement a
    * serialize() function or not, as well as objects that are iterable.
    * */
  template<class Archive, class T, bool HasSerializer, bool IsIterable>
  struct wrapper;
}

/** is_iterable<T>::value is true if T is iterable, false otherwise */
template <typename T>
using is_iterable = decltype(detail::is_iterable_impl<T>(0));

/** is_serializable<T>::value is true if T declares serializable, false
 * otherwise */
template <typename T>
using is_serializable = decltype(detail::is_serializable_impl<T>(0));

/** is_true_iterable<T>: iterable and NOT string */
template<typename T>
struct is_true_iterable
{
  const static bool value = is_iterable<T>::value &&
                            !std::is_same<T, std::string>::value;
};

// =============================================================================
// Traits for fundamental types

/** Type traits to differentiate iterable types */
template<class T, bool IsIterable> struct type_name_impl
{};

/** Small type traits class
  *
  * We use this small struct to generate strings corresponding to types at
  * compile time, e.g. double -> "double", in a somewhat platform independent
  * way. If you get error involving this class this probably means that you need
  * to register a new type: add the corresponding REGISTER_TYPE_NAME or
  * REGISTER_TYPE_NAME_OTHER instruction in serialization.cpp or below your
  * class definition.
  * */
template<class T>
struct type_name
  : type_name_impl<T, is_true_iterable<T>::value>
{};

/** Iterable types are all called 'array(...)'
  *
  * Iterable types must have T::value_type, as in the standard containers, to
  * declare what type they are storing. Remember that in C++ the types are
  * homogeneous within a container.
  * */
template<class T>
struct type_name_impl<T, true>
{
  static std::string name()
  { return "array("+type_name<typename T::value_type>::name()+")"; }
};

/** Register type X to have name Y */
#define REGISTER_TYPE_NAME_OTHER(X, Y) \
  template<> struct type_name_impl<X, false> \
  { static std::string name() { return Y; } };
/** Register type X to have name "X" */
#define REGISTER_TYPE_NAME(X) REGISTER_TYPE_NAME_OTHER(X, #X)

// We register the names of the different fundamental types
REGISTER_TYPE_NAME(bool)
REGISTER_TYPE_NAME(int)
REGISTER_TYPE_NAME(unsigned)
REGISTER_TYPE_NAME(short)
REGISTER_TYPE_NAME(unsigned short)
REGISTER_TYPE_NAME(long)
REGISTER_TYPE_NAME(unsigned long)
REGISTER_TYPE_NAME(long long)
REGISTER_TYPE_NAME(unsigned long long)
REGISTER_TYPE_NAME(float)
REGISTER_TYPE_NAME(double)
REGISTER_TYPE_NAME(long double)
REGISTER_TYPE_NAME_OTHER(std::string, "string")
REGISTER_TYPE_NAME_OTHER(char, "string")
REGISTER_TYPE_NAME_OTHER(wchar_t, "string")
REGISTER_TYPE_NAME_OTHER(char16_t, "string")
REGISTER_TYPE_NAME_OTHER(char32_t, "string")
REGISTER_TYPE_NAME_OTHER(const char*, "string")
REGISTER_TYPE_NAME_OTHER(const unsigned char*, "string")
REGISTER_TYPE_NAME_OTHER(const wchar_t*, "string")
REGISTER_TYPE_NAME_OTHER(const char16_t*, "string")
REGISTER_TYPE_NAME_OTHER(const char32_t*, "string")

// =============================================================================
// Traits to obtain shape of iterables

/** Type traits to differentiate iterable types */
template<class T, bool IsIterable> struct get_shape_impl
{};

template<class T>
std::vector<size_t> get_shape(const T& iterable)
{
  auto ret = get_shape_impl<T, is_true_iterable<T>::value>::shape(iterable);
  std::reverse(ret.begin(), ret.end());
  return ret;
}

template<class T>
struct get_shape_impl<T, false>
{
  static std::vector<size_t> shape(const T& iterable)
  {
    return {};
  }
};

template<class T>
struct get_shape_impl<T, true>
{
  static std::vector<size_t> shape(const T& iterable)
  {
    using U = typename std::decay<decltype(*std::begin(iterable))>::type;
    auto ret = get_shape_impl<
                 U, is_true_iterable<U>::value
                 >::shape(*std::begin(iterable));
    ret.push_back(iterable.size());
    return ret;
  }
};

// =============================================================================
// The archive types

/** Output json archive */
class oarchive : public std::ostream
{
  /** The output stream */
  std::ostream& stream;
  /** Indentation level (we aren't monsters) */
  unsigned level = 0;
  /** Add indentation level */
  void indent(unsigned=1);
  /** Remove indentation level */
  void unindent(unsigned=1);
  /** Used to avoid having commas everywhere */
  bool first = true;

  /** Open new group { ... } */
  void open_group(const std::string& = "{");
  /** Close current group */
  void close_group(const std::string& = "}");
  /** Create a new line, taking care of the commas and indentation */
  void new_line();

  /** Write key ( "key" : element ) */
  void add_key(const std::string&);
  /** Write element ( "key" : element ), while putting quotes for strings */
  template<class T>
  void add_element(const T&);
  /** Write key/element pair */
  template<class T>
  void add(const std::string& key, const T& element)
  { add_key(key); add_element(element); }
  /** Write iterable element -> [ ... ] */
  template<class T>
  void add_iterable(const T&);

  /** Bad value flag (nans) */
  bool f_bad_value = false;
public:
  /** Contructor
   *
   * Arguments are the output stream, the id, and version number.
   * */
  oarchive(std::ostream&, std::string, unsigned=0);
  /** Write final characters to stream and destruct */
  ~oarchive();

  /** The archive output operator */
  template<class T>
  oarchive& operator&(const std::pair<T&, std::string>&);

  template<class Archive, class T, bool HasSerializer, bool IsIterable>
  friend struct detail::wrapper;

  template<class T>
  void serialize(const T&);

  /** Return true if a nan was found while writting */
  bool bad_value() const { return f_bad_value; }
};

/** For non-string types: put into stream using stringstream directly */
template<class T>
void oarchive::add_element(const T& value)
{
  // detect nans
  if(value!=value) f_bad_value = true;
  // write
  std::stringstream ss;
  ss << value;
  stream << ss.str();
}
/** Template specialization for const char*: add enclsoing braces */
template<>
void oarchive::add_element<const char*>(const char* const&);
/** Template specialization for strings: add enclosing braces */
template<>
void oarchive::add_element<std::string>(const std::string&);

template<class T>
void oarchive::add_iterable(const T& value)
{
  stream << "[ ";
  bool f = true;
  for(const auto& i : value)
  {
    // commas
    if(f) f=false;
    else stream << ", ";
    // add elements
    serialize(i);
  }
  stream << " ]";
}

/** Input archive
  *
  * Not implemented.
  */
struct iarchive
{
  /** The archive input operator */
  template<class T>
  iarchive& operator&(const std::pair<T&, std::string>&)
  { return *this; }
};

// =============================================================================
// Template specializations for traits

namespace detail
{
  /** If object has a serialize function */
  template<class Archive, class T, bool IsIterable>
  struct wrapper<Archive, T, true, IsIterable>
  {
    wrapper(Archive& ar, const T& obj)
    {
      ar.open_group();
      // this cast is fine because we are only writting
      const_cast<T&>(obj).serialize(ar);
      ar.close_group();
    }
  };

  /** If the object has no serialize function (non-iterable) */
  template<class Archive, class T>
  struct wrapper<Archive, T, false, false>
  {
    wrapper(Archive& ar, const T& obj)
    {
      ar.add_element(obj);
    }
  };

  /** If the object has no serialize function (iterable) */
  template<class Archive, class T>
  struct wrapper<Archive, T, false, true>
  {
    wrapper(Archive& ar, const T& obj)
    {
      ar.add_iterable(obj);
    }
  };
}

// =============================================================================
// Here is where things happen

/** Most of the magic is happening here
  *
  * This function is merely a shortcut for calling the wrapper.
  * */
template<class T>
void oarchive::serialize(const T& t)
{
  detail::wrapper
    < oarchive, T,
      is_serializable<T>::value,
      is_true_iterable<T>::value
    > (*this, t);
}

template<class T>
oarchive& oarchive::operator&(const std::pair<T&, std::string>& t)
{
  // open a new group with variable name
  add_key(t.second);
  open_group();
  // write the type
  add_key("type");
  add_element(type_name<T>::name());
  // in the case it is iterable, print shape
  if(is_true_iterable<T>::value)
  {
    add_key("shape");
    serialize(get_shape(t.first));
  }
  // write value
  add_key("value");
  // serialize the object
  serialize(t.first);
  // close the group
  close_group();
  // bye
  return *this;
}

#endif//SERIALIZATION_HPP_
