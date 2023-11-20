/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_EMBED_INDENT_H
#define INCLUDE_SPARROW_EMBED_INDENT_H

#include <string>

template<typename T>
struct Indented {
  constexpr static unsigned width = 2;

  unsigned startLevel;
  const T& bound;

  std::string operator()(const unsigned level) const {
    return std::string(width * (startLevel + level), ' ');
  }

  Indented<T> increment() const {
    return Indented<T>{startLevel + 1, bound};
  }
};

struct Indent {
  constexpr static unsigned width = 2;

  static std::string level(const unsigned a) {
    return std::string(width * a, ' ');
  }

  template<typename T>
  static auto bind(unsigned level, const T& t) {
    return Indented<T>{level, t};
  }
};

#endif
