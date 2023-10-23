//
// Created by Ludwig Schmid on 19.10.23.
//

#pragma once

#include "namap/NeutralAtomArchitecture.hpp"
#include "namap/NeutralAtomDefinitions.hpp"
#include "namap/NeutralAtomUtils.hpp"

namespace qc {

// Class to manage hardware qubit handling
class HardwareQubits {
protected:
  //        std::map<Qubit, Qubit>      circToHw;
  std::map<HwQubit, CoordIndex> hwToCoordIdx;
  SymmetricMatrix               swapDistances;

  void initSwapDistances(const NeutralAtomArchitecture& arch);
  void updateSwapDistances(const NeutralAtomArchitecture& arch, HwQubit qubit);

public:
  HardwareQubits() = delete;
  HardwareQubits(const NeutralAtomArchitecture& arch,
                 InitialCoordinateMapping&      initialCoordinateMapping)
      : swapDistances(arch.getNqubits()) {
    switch (initialCoordinateMapping) {
    case Trivial:
      for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
        hwToCoordIdx.insert({i, i});
      }
      break;
    }
    initSwapDistances(arch);
  }

  [[nodiscard]] inline fp getSwapDistance(HwQubit q1, HwQubit q2) const {
    return swapDistances(q1, q2);
  }
  [[nodiscard]] inline CoordIndex getCoordIndex(HwQubit qubit) const {
    return hwToCoordIdx.at(qubit);
  }
  [[nodiscard]] inline Qubit getQubit(CoordIndex coordIndex) const {
    for (auto const& [qubit, index] : hwToCoordIdx) {
      if (index == coordIndex) {
        return qubit;
      }
    }
    throw std::runtime_error("There is no qubit at this coordinate " +
                             std::to_string(coordIndex));
  }

  void move(HwQubit hwQubit, CoordIndex newCoord,
            NeutralAtomArchitecture& arch);

  std::vector<Swap>    getNearbySwaps(HwQubit q);
  std::vector<HwQubit> getNearbyQubits(HwQubit q);
  fp                   getTotalDistance(std::set<HwQubit>& qubits);
  std::set<HwQubit>    getBlockedQubits(const std::set<HwQubit>&       qubits,
                                        const NeutralAtomArchitecture& arch);
};
} // namespace qc
