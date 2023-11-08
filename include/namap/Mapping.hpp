//
// Created by Ludwig Schmid on 19.10.23.
//

#pragma once

#include "namap/NeutralAtomDefinitions.hpp"
#include "namap/NeutralAtomUtils.hpp"

namespace qc {

// class to manage the mapping between circuit qubits and hardware qubits
// in a bijective manner
class Mapping {
protected:
  // std::map<Qubit, HwQubit>
  Permutation circToHw;

public:
  Mapping() = default;
  Mapping(size_t nQubits, InitialMapping initialMapping) {
    switch (initialMapping) {
    case Identity:
      for (size_t i = 0; i < nQubits; ++i) {
        circToHw.insert({i, i});
      }
      break;
    }
  }
  void inline setCircuitQubit(Qubit qubit, HwQubit hwQubit) {
    circToHw[qubit] = hwQubit;
  }

  [[nodiscard]] inline HwQubit getHwQubit(Qubit qubit) const {
    return circToHw.at(qubit);
  }

  [[nodiscard]] inline std::set<HwQubit>
  getHwQubits(std::set<Qubit>& qubits) const {
    std::set<HwQubit> hwQubits;
    for (const auto& qubit : qubits) {
      hwQubits.insert(this->getHwQubit(qubit));
    }
    return hwQubits;
  }

  [[nodiscard]] inline Qubit getCircQubit(HwQubit qubit) const {
    for (const auto& [circQubit, hwQubit] : circToHw) {
      if (hwQubit == qubit) {
        return circQubit;
      }
    }
    throw std::runtime_error("Hardware qubit: " + std::to_string(qubit) +
                             " not found in mapping");
  }

  [[nodiscard]] inline bool isMapped(HwQubit qubit) const {
    return std::any_of(
        circToHw.begin(), circToHw.end(),
        [qubit](const auto& pair) { return pair.second == qubit; });
  }

  inline void mapToHwQubits(Operation* op) const {
    op->setTargets(circToHw.apply(op->getTargets()));
    op->setControls(circToHw.apply(op->getControls()));
  }

  void swap(Swap swap);
};

} // namespace qc
