//
// Created by Ludwig Schmid on 19.10.23.
//

#include "namap/HardwareQubits.hpp"

#include "namap/Mapping.hpp"

namespace qc {
void HardwareQubits::initSwapDistances(const NeutralAtomArchitecture& arch) {
  swapDistances = SymmetricMatrix(arch.getNqubits());
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      swapDistances(i, j) =
          arch.getSwapDistance(hwToCoordIdx.at(i), hwToCoordIdx.at(j));
    }
  }
}

void HardwareQubits::updateSwapDistances(const NeutralAtomArchitecture& arch,
                                         HwQubit                        qubit) {
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    swapDistances(i, qubit) =
        arch.getSwapDistance(hwToCoordIdx.at(i), hwToCoordIdx.at(qubit));
  }
}

void HardwareQubits::move(HwQubit hwQubit, CoordIndex newCoord,
                          NeutralAtomArchitecture& arch) {
  if (newCoord >= arch.getNpositions()) {
    throw std::runtime_error("Invalid coordinate");
  }
  // check if new coordinate is already occupied
  for (const auto& [qubit, coord] : hwToCoordIdx) {
    if (coord == newCoord) {
      throw std::runtime_error("Coordinate already occupied");
    }
  }
  hwToCoordIdx.at(hwQubit) = newCoord;
  updateSwapDistances(arch, hwQubit);
}

std::vector<Swap> HardwareQubits::getNearbySwaps(qc::HwQubit q) {
  std::vector<Swap> swaps;
  auto              nearbyQubits = getNearbyQubits(q);
  swaps.reserve(nearbyQubits.size());
  for (const auto& nearbyQubit : nearbyQubits) {
    swaps.emplace_back(q, nearbyQubit);
  }
  return swaps;
}

std::vector<HwQubit> HardwareQubits::getNearbyQubits(qc::HwQubit q) {
  std::vector<HwQubit> nearbyQubits;
  for (uint32_t i = 0; i < swapDistances.getSize(); ++i) {
    if (i == q) {
      continue;
    }
    if (swapDistances(q, i) == 0) {
      nearbyQubits.emplace_back(i);
    }
  }
  return nearbyQubits;
}

fp HardwareQubits::getTotalDistance(std::set<Qubit>& qubits) {
  // two qubit gates
  if (qubits.size() == 2) {
    auto it = qubits.begin();
    auto q1 = *it;
    auto q2 = *(++it);
    return swapDistances(q1, q2);
  }
  if (qubits.size() == 3) {
    // TODO substitute with special case taking into consideration the geometry
    auto it = qubits.begin();
    auto q1 = *it;
    auto q2 = *(++it);
    auto q3 = *(++it);
    return swapDistances(q1, q2) + swapDistances(q2, q3) +
           swapDistances(q1, q3);
  }
  // more than three qubits just minimize total distance
  fp totalDistance = 0;
  for (auto it1 = qubits.begin(); it1 != qubits.end(); ++it1) {
    for (auto it2 = std::next(it1); it2 != qubits.end(); ++it2) {
      totalDistance += swapDistances(*it1, *it2);
    }
  }
  return totalDistance;
}
} // namespace qc
