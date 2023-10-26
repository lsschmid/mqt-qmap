//
// Created by Ludwig Schmid on 19.10.23.
//

#include "namap/HardwareQubits.hpp"

#include "namap/Mapping.hpp"

namespace qc {
void HardwareQubits::initCompactSwapDistances(
    const NeutralAtomArchitecture& arch) {
  // only valid for trivial/compact initial layout
  swapDistances = SymmetricMatrix(arch.getNqubits());
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      swapDistances(i, j) =
          arch.getSwapDistance(hwToCoordIdx.at(i), hwToCoordIdx.at(j));
    }
  }
}

void HardwareQubits::initNearbyQubits(const qc::NeutralAtomArchitecture& arch) {
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    computeNearbyQubits(i);
  }
}

void HardwareQubits::computeSwapDistance(HwQubit q1, HwQubit q2) {
  std::queue<HwQubit>  q;
  std::vector<bool>    visited(swapDistances.getSize(), false);
  std::vector<HwQubit> parent(swapDistances.getSize(), q2);

  q.push(q1);
  visited[q1] = true;
  parent[q1]  = q1;
  bool found  = false;
  while (!q.empty() && !found) {
    auto current = q.front();
    q.pop();
    for (const auto& nearbyQubit : nearbyQubits.at(current)) {
      if (!visited[nearbyQubit]) {
        q.push(nearbyQubit);
        visited[nearbyQubit] = true;
        parent[nearbyQubit]  = current;
        if (nearbyQubit == q2) {
          found = true;
          break;
        }
      }
    }
  }
  if (!found) {
    swapDistances(q1, q2) = std::numeric_limits<fp>::infinity();
    return;
  }
  // recreate path
  std::vector<HwQubit> path;
  auto                 current = q2;
  while (current != q1) {
    path.emplace_back(current);
    current = parent[current];
  }
  path.emplace_back(q1);
  // update swap distances along path
  for (uint32_t start = 0; start < path.size() - 1; ++start) {
    for (uint32_t end = start + 1; end < path.size(); ++end) {
      swapDistances(path[start], path[end]) = end - start - 1;
    }
  }
  //
  //  swapDistances(q1, q2) = static_cast<fp>(path.size() - 2);
}

void HardwareQubits::updateSwapDistances(
    const qc::NeutralAtomArchitecture& arch, qc::HwQubit qubit) {
  // reset swap distances
  swapDistances = SymmetricMatrix(arch.getNqubits(), -1);
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

  // remove qubit from old nearby qubits
  auto prevNearbyQubits = nearbyQubits.at(hwQubit);
  for (const auto& qubit : prevNearbyQubits) {
    nearbyQubits.at(qubit).erase(std::find(
        nearbyQubits.at(qubit).begin(), nearbyQubits.at(qubit).end(), hwQubit));
  }
  // move qubit and compute new nearby qubits
  hwToCoordIdx.at(hwQubit) = newCoord;
  computeNearbyQubits(hwQubit);

  // add qubit to new nearby qubits
  auto newNearbyQubits = nearbyQubits.at(hwQubit);
  for (const auto& qubit : newNearbyQubits) {
    nearbyQubits.at(qubit).emplace_back(hwQubit);
  }

  // update/reset swap distances
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

void HardwareQubits::computeNearbyQubits(qc::HwQubit q) {
  std::vector<HwQubit> newNearbyQubits;
  auto                 coordQ = hwToCoordIdx.at(q);
  for (const auto& coord : hwToCoordIdx) {
    if (coord.first == q) {
      continue;
    }
    if (arch.getEuclidianDistance(coordQ, coord.second) <=
        arch.getInteractionRadius()) {
      newNearbyQubits.emplace_back(coord.first);
    }
  }
  nearbyQubits.insert_or_assign(q, newNearbyQubits);
}

fp HardwareQubits::getTotalDistance(std::set<HwQubit>& hwQubits) {
  // two qubit gates
  if (hwQubits.size() == 2) {
    auto it = hwQubits.begin();
    auto q1 = *it;
    auto q2 = *(++it);
    return getSwapDistance(q1, q2);
  }
  if (hwQubits.size() == 3) {
    // TODO substitute with special case taking into consideration the
    // geometry
    auto it = hwQubits.begin();
    auto q1 = *it;
    auto q2 = *(++it);
    auto q3 = *(++it);
    return getSwapDistance(q1, q2) + getSwapDistance(q2, q3) +
           getSwapDistance(q1, q3);
  }
  // more than three hwQubits just minimize total distance
  fp totalDistance = 0;
  for (auto it1 = hwQubits.begin(); it1 != hwQubits.end(); ++it1) {
    for (auto it2 = std::next(it1); it2 != hwQubits.end(); ++it2) {
      totalDistance += getSwapDistance(*it1, *it2);
    }
  }
  return totalDistance;
}

std::set<HwQubit>
HardwareQubits::getBlockedQubits(const std::set<HwQubit>&           qubits,
                                 const qc::NeutralAtomArchitecture& arch) {
  std::set<HwQubit> blockedQubits;
  for (const auto& qubit : qubits) {
    for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
      if (i == qubit) {
        continue;
      }
      // do a preselection
      // now check exact difference
      auto const distance =
          arch.getEuclidianDistance(hwToCoordIdx.at(qubit), hwToCoordIdx.at(i));
      if (distance <= arch.getBlockingFactor() * arch.getInteractionRadius()) {
        blockedQubits.insert(i);
      }
    }
  }
  return blockedQubits;
}

} // namespace qc
