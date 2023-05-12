/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_EMBED_COMMON_H
#define INCLUDE_SPARROW_EMBED_COMMON_H

#include "boost/optional.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

inline std::ostream& sep(std::ostream& os) {
  os << ", ";
  return os;
}

inline std::ostream& nl(std::ostream& os) {
  os << "\n";
  return os;
}

inline std::string sp(unsigned i) {
  return std::string(i, ' ');
}

struct CppFileStructure {
  std::vector<std::string> includes;
  std::vector<std::string> namespaces;
  boost::optional<std::string> guard = boost::none;

  void writeHeader(std::ostream& os) {
    /* We need to obfuscate the doxygen copyright statement from our license
     * checker script, which allows only a single such statement in each file.
     */
    os << R"delim(/**
 * @file
 * @copy)delim"
       << R"delim(right This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 *
 * @note This file was generated by the Embed binary from runtime values.
 *   Prefer improving the generator over editing this file whenever possible.
 *
)delim";

    // Write some more details for cpp files
    if (!guard) {
      os << R"delim( * This file contains functions generating runtime values. It was
 * generated from runtime values of its return type. It is not intended to be
 * human-readable. A small guide: Return values are directly brace-initialized
 * in deep nesting to keep file size to a minimum. Types are annotated only when
 * necessary. Floating point values are represented in hexadecimal (see
 * std::hexfloat) to ensure serialization does not cause loss of accuracy.
 *
 * The functions defined here might be declared and called elsewhere entirely.
 *
)delim";
    }

    os << R"delim( */
)delim";

    if (guard) {
      os << "#ifndef " << guard.value() << nl;
      os << "#define " << guard.value() << nl << nl;
    }

    for (auto& include : includes) {
      if (include[0] == '<' || include[0] == '\"') {
        os << "#include " << include << nl;
      }
      else {
        os << "#include \"" << include << "\"" << nl;
      }
    }
    os << nl;

    for (auto& name : namespaces) {
      os << "namespace " << name << " {" << nl;
    }

    os << nl;
  }

  void writeFooter(std::ostream& os) {
    os << nl;
    for (auto iter = std::rbegin(namespaces); iter != std::rend(namespaces); ++iter) {
      os << "} // namespace " << *iter << nl;
    }

    if (guard) {
      os << nl << "#endif";
    }
  }
};

inline std::string lower(std::string a) {
  std::transform(std::begin(a), std::end(a), std::begin(a), [](unsigned char c) { return std::tolower(c); });
  return a;
}

inline std::string upper(std::string a) {
  std::transform(std::begin(a), std::end(a), std::begin(a), [](unsigned char c) { return std::toupper(c); });
  return a;
}

namespace detail {

template<typename Arg>
void writePack(std::ostream& os, Arg arg) {
  os << arg;
}

template<typename Arg, typename... Args>
void writePack(std::ostream& os, Arg arg1, Args... args) {
  os << arg1 << sep;
  writePack(os, args...);
}

template<typename Tuple, std::size_t... Inds>
void writeTuplePackHelper(std::ostream& os, const Tuple& tup, std::index_sequence<Inds...> /* inds */
) {
  writePack(os, std::get<Inds>(tup)...);
}

template<typename Tuple>
void writeTuplePack(std::ostream& os, const Tuple& tup) {
  writeTuplePackHelper(os, tup, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
}

} // namespace detail

template<typename... Args>
struct Joiner {
  Joiner(const Args&... args) : bound(args...) {
  }

  std::tuple<const Args&...> bound;
};

template<typename... Args>
auto join(const Args&... args) {
  return Joiner<Args...>(args...);
}

template<typename... Args>
std::ostream& operator<<(std::ostream& os, const Joiner<Args...>& join) {
  detail::writeTuplePack(os, join.bound);
  return os;
}

/* Note: We need to wrap optionals because boost provides its own output
 * operators that don't serve our purpose. We can't use the type directly
 * because there's already a declaration available for it by inclusion of boost
 * optional.
 */
template<typename T>
struct OptionalWrapper {
  boost::optional<T> optional;
};

template<typename T>
auto wrapOptional(const boost::optional<T>& optional) {
  return OptionalWrapper<T>{optional};
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const OptionalWrapper<T>& o) {
  if (o.optional) {
    os << o.optional.value();
  }
  else {
    os << "boost::none";
  }
  return os;
}

#endif
