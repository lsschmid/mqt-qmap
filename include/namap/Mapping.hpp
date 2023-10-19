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
  std::map<Qubit, HwQubit> circToHw;
  std::map<HwQubit, Qubit> hwToCirc;

public:
  Mapping() = default;
  Mapping(size_t nQubits, InitialMapping initialMapping) {
    switch (initialMapping) {
    case Identity:
      for (size_t i = 0; i < nQubits; ++i) {
        circToHw.insert({i, i});
        hwToCirc.insert({i, i});
      }
      break;
    }
  }
  void inline setCircuitQubit(Qubit qubit, HwQubit hwQubit) {
    circToHw.insert({qubit, hwQubit});
    hwToCirc.insert({hwQubit, qubit});
  }
  void inline removeCircuitQubit(Qubit qubit) {
    auto hwQubit = circToHw.at(qubit);
    circToHw.erase(qubit);
    hwToCirc.erase(hwQubit);
  }

  [[nodiscard]] inline HwQubit getHwQubit(Qubit qubit) const {
    return circToHw.at(qubit);
  }
  [[nodiscard]] inline Qubit getCircQubit(HwQubit qubit) const {
    return hwToCirc.at(qubit);
  }
  [[nodiscard]] inline bool isMapped(HwQubit qubit) const {
    return hwToCirc.find(qubit) != hwToCirc.end();
  }

  void swap(Swap swap);
};

} // namespace qc
