/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */
#ifndef INCLUDE_SPARROW_CONSTEXPR_NDDO_PARAMETERS_H
#define INCLUDE_SPARROW_CONSTEXPR_NDDO_PARAMETERS_H

#include "Sparrow/Implementations/Nddo/Parameters.h"
#include <algorithm>
#include <array>
#include <tuple>

namespace Scine {
namespace Sparrow {
namespace nddo {

template<std::size_t s>
struct ConstexprAtomic {
  static constexpr std::size_t size = s;

  int Z;
  Parameters::Atomic::Pack pack;
  std::array<Parameters::Atomic::GaussianRepulsion, size> gaussianRepulsion;
};

template<typename... GRep>
constexpr auto atomic(int Z, Parameters::Atomic::Pack pack, GRep&&... greps) {
  return ConstexprAtomic<sizeof...(greps)>{Z, std::move(pack), {std::forward<GRep>(greps)...}};
}

template<std::size_t s>
Parameters::Atomic runtime(const ConstexprAtomic<s>& a) {
  Parameters::Atomic atomic;
  atomic.pack = a.pack;
  atomic.gaussianRepulsion.resize(s);
  std::copy(std::begin(a.gaussianRepulsion), std::end(a.gaussianRepulsion), std::begin(atomic.gaussianRepulsion));
  return atomic;
}

namespace detail {

template<typename AtomicsTuple, std::size_t... Inds>
std::unordered_map<int, Parameters::Atomic> runtime(AtomicsTuple&& tup, std::index_sequence<Inds...> /* inds */) {
  return std::unordered_map<int, Parameters::Atomic>{
      {std::get<Inds>(tup).Z, runtime(std::get<Inds>(tup))}...,
  };
}

} // namespace detail

template<typename AtomicsTuple>
std::unordered_map<int, Parameters::Atomic> runtime(AtomicsTuple&& tup) {
  return detail::runtime(std::forward<AtomicsTuple>(tup),
                         std::make_index_sequence<std::tuple_size<std::decay_t<AtomicsTuple>>::value>{});
}

struct ConstexprDiatomic {
  int Z1 = 0;
  int Z2 = 0;
  Parameters::Diatomic constants;
};

constexpr ConstexprDiatomic diatomic(int Z1, int Z2, double exponent, double factor) {
  return {Z1, Z2, {exponent, factor}};
}

template<typename AtomicList, std::size_t d>
struct ConstexprParameters {
  static constexpr std::size_t diatomicSize = d;

  AtomicList atomic;
  std::array<ConstexprDiatomic, d> diatomic;
};

template<typename... Atomics>
constexpr auto collect(Atomics&&... atomics) {
  return std::make_tuple(std::forward<Atomics>(atomics)...);
}

template<typename AtomicList, typename... Diatomics>
constexpr auto make_parameters(AtomicList&& atomics, Diatomics... diatomics) {
  return ConstexprParameters<AtomicList, sizeof...(Diatomics)>{std::forward<AtomicList>(atomics),
                                                               {
                                                                   std::forward<Diatomics>(diatomics)...,
                                                               }};
}

template<typename AtomicList, std::size_t s>
Parameters runtime(const ConstexprParameters<AtomicList, s>& a) {
  Parameters nddo;

  nddo.atomic = runtime(a.atomic);

  for (const ConstexprDiatomic& diatomic : a.diatomic) {
    nddo.diatomic.emplace(std::make_pair(diatomic.Z1, diatomic.Z2), diatomic.constants);
  }

  return nddo;
}

} // namespace nddo
} // namespace Sparrow
} // namespace Scine

#endif
