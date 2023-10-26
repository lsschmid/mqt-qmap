//
// Created by Ludwig Schmid on 19.10.23.
//

#pragma once

#include "limits"
#include "namap/NeutralAtomArchitecture.hpp"
#include "namap/NeutralAtomDefinitions.hpp"
#include "namap/NeutralAtomUtils.hpp"
#include "random"

namespace qc {

// Class to manage hardware qubit handling
class HardwareQubits {
protected:
  //        std::map<Qubit, Qubit>      circToHw;
  NeutralAtomArchitecture                 arch;
  std::map<HwQubit, CoordIndex>           hwToCoordIdx;
  SymmetricMatrix                         swapDistances;
  std::map<HwQubit, std::vector<HwQubit>> nearbyQubits;

  void initCompactSwapDistances(const NeutralAtomArchitecture& arch);
  void initNearbyQubits(const NeutralAtomArchitecture& arch);
  void updateSwapDistances(const NeutralAtomArchitecture& arch, HwQubit qubit);
  void computeSwapDistance(HwQubit q1, HwQubit q2);
  void computeNearbyQubits(HwQubit qubit);

public:
  HardwareQubits() = delete;
  HardwareQubits(const NeutralAtomArchitecture& arch,
                 InitialCoordinateMapping&      initialCoordinateMapping)
      : arch(arch), swapDistances(arch.getNqubits()) {
    switch (initialCoordinateMapping) {
    case Trivial:
      for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
        hwToCoordIdx.insert({i, i});
      }
      initCompactSwapDistances(arch);
      break;
    case Random:
      std::vector<CoordIndex> indices(arch.getNpositions());
      std::iota(indices.begin(), indices.end(), 0);
      std::random_device rd;
      std::mt19937       g(rd());
      std::shuffle(indices.begin(), indices.end(), g);
      for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
        hwToCoordIdx.insert({i, indices[i]});
      }

      swapDistances = SymmetricMatrix(arch.getNqubits(), -1);
    }
    initNearbyQubits(arch);
  }

  [[nodiscard]] inline fp getSwapDistance(HwQubit q1, HwQubit q2) {
    if (swapDistances(q1, q2) < 0) {
      computeSwapDistance(q1, q2);
    }
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
  [[nodiscard]] inline std::vector<HwQubit> getNearbyQubits(HwQubit q) const {
    return nearbyQubits.at(q);
  }

  void move(HwQubit hwQubit, CoordIndex newCoord,
            NeutralAtomArchitecture& arch);

  std::vector<Swap> getNearbySwaps(HwQubit q);
  fp                getTotalDistance(std::set<HwQubit>& qubits);
  std::set<HwQubit> getBlockedQubits(const std::set<HwQubit>&       qubits,
                                     const NeutralAtomArchitecture& arch);
};
} // namespace qc
