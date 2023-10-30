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

  void initCompactSwapDistances();
  void initNearbyQubits();
  void updateSwapDistances(HwQubit qubit);
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
      initCompactSwapDistances();
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
    initNearbyQubits();
  }

  [[nodiscard]] inline fp getSwapDistance(HwQubit q1, HwQubit q2,
                                          bool closeBy = true) {
    if (q1 == q2) {
      return 0;
    }
    if (swapDistances(q1, q2) < 0) {
      computeSwapDistance(q1, q2);
    }
    if (closeBy) {
      return swapDistances(q1, q2);
    }
    return swapDistances(q1, q2) + 1;
  }
  [[nodiscard]] fp getSwapDistanceMove(CoordIndex idx, HwQubit target);
  [[nodiscard]] inline CoordIndex getCoordIndex(HwQubit qubit) const {
    return hwToCoordIdx.at(qubit);
  }
  [[nodiscard]] inline std::set<CoordIndex>
  getCoordIndices(std::set<HwQubit>& hwQubits) const {
    std::set<CoordIndex> coordIndices;
    for (auto const& hwQubit : hwQubits) {
      coordIndices.insert(this->getCoordIndex(hwQubit));
    }
    return coordIndices;
  }
  [[nodiscard]] inline HwQubit getHwQubit(CoordIndex coordIndex) const {
    for (auto const& [hwQubit, index] : hwToCoordIdx) {
      if (index == coordIndex) {
        return hwQubit;
      }
    }
    throw std::runtime_error("There is no qubit at this coordinate " +
                             std::to_string(coordIndex));
  }
  [[nodiscard]] inline std::vector<HwQubit> getNearbyQubits(HwQubit q) const {
    return nearbyQubits.at(q);
  }
  [[nodiscard]] inline bool isMapped(CoordIndex idx) const {
    return !std::none_of(
        hwToCoordIdx.begin(), hwToCoordIdx.end(),
        [idx](const auto& pair) { return pair.second == idx; });
  }
  [[nodiscard]] inline std::set<CoordIndex>
  getNearbyCoordinates(HwQubit q) const {
    return this->arch.getNearbyCoordinates(this->getCoordIndex(q));
  }

  void move(HwQubit hwQubit, CoordIndex newCoord);

  std::vector<Swap>       getNearbySwaps(HwQubit q);
  std::set<CoordIndex>    getNearbyFreeCoordinates(HwQubit q);
  fp                      getTotalDistance(std::set<HwQubit>& qubits);
  std::set<HwQubit>       getBlockedQubits(const std::set<HwQubit>& qubits);
  std::vector<CoordIndex> findClosestFreeCoord(HwQubit               qubit,
                                               MoveVector::Direction direction);
};
} // namespace qc
