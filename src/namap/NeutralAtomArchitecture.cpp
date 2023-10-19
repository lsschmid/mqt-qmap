//
// Created by Ludwig Schmid on 13.10.23.
//

#include "namap/NeutralAtomArchitecture.hpp"

#include "nlohmann/json.hpp"
#include "operations/OpType.hpp"

namespace qc {

void NeutralAtomArchitecture::loadJson(const std::string& filename) {
  nlohmann::json jsonData;
  std::ifstream  architecureFile(filename);

  if (!architecureFile.is_open()) {
    throw std::runtime_error("Could not open file " + filename);
  }
  try {
    architecureFile >> jsonData;
    architecureFile.close();

    // Load properties
    nlohmann::json jsonDataProperties = jsonData["properties"];
    this->properties                  = Properties(
        jsonDataProperties["nRows"], jsonDataProperties["nColumns"],
        jsonDataProperties["nAods"], jsonDataProperties["nAodCoordinates"],
        jsonDataProperties["interQubitDistance"],
        jsonDataProperties["interactionRadius"],
        jsonDataProperties["blockingFactor"]);

    // Load parameters
    const nlohmann::json jsonDataParameters = jsonData["parameters"];
    this->parameters                        = Parameters();
    this->parameters.nQubits                = jsonDataParameters["nQubits"];

    std::map<OpType, fp> gateTimes;
    for (const auto& [key, value] : jsonDataParameters["gateTimes"].items()) {
      gateTimes.insert({OP_NAME_TO_TYPE.at(key), value});
    }
    this->parameters.gateTimes = gateTimes;
    std::map<OpType, fp> gateAverageFidelities;
    for (const auto& [key, value] :
         jsonDataParameters["gateAverageFidelities"].items()) {
      gateAverageFidelities.insert({OP_NAME_TO_TYPE.at(key), value});
    }
    this->parameters.gateAverageFidelities = gateAverageFidelities;
    std::map<ShuttlingType, fp> shuttlingTimes;

    for (const auto& [key, value] :
         jsonDataParameters["shuttlingTimes"].items()) {
      shuttlingTimes.insert({SHUTTLING_OP_NAME_TO_TYPE.at(key), value});
    }
    this->parameters.shuttlingTimes = shuttlingTimes;
    std::map<ShuttlingType, fp> shuttlingAverageFidelities;
    for (const auto& [key, value] :
         jsonDataParameters["shuttlingAverageFidelities"].items()) {
      shuttlingAverageFidelities.insert(
          {SHUTTLING_OP_NAME_TO_TYPE.at(key), value});
    }
    this->parameters.shuttlingAverageFidelities = shuttlingAverageFidelities;

    this->parameters.decoherenceTimes =
        NeutralAtomArchitecture::Parameters::DecoherenceTimes(
            jsonDataParameters["decoherenceTimes"]["t1"],
            jsonDataParameters["decoherenceTimes"]["t2"]);

  } catch (std::exception& e) {
    throw std::runtime_error("Could not parse JSON file " + filename + ": " +
                             e.what());
  }

  // apply changes to the object
  this->name = jsonData["name"];

  this->createCoordinates();
  this->computeSwapDistances(this->properties.getInteractionRadius());
}
void NeutralAtomArchitecture::createCoordinates() {
  for (std::uint16_t i = 0; i < this->properties.getNpositions(); i++) {
    this->coordinates.emplace_back(i % this->properties.getNcolumns(),
                                   i / this->properties.getNcolumns());
  }
}
NeutralAtomArchitecture::NeutralAtomArchitecture(const std::string& filename) {
  this->loadJson(filename);
}

void NeutralAtomArchitecture::computeSwapDistances(fp interactionRadius) {
  // compute diagonal distances
  struct DiagonalDistance {
    std::uint32_t x;
    std::uint32_t y;
    fp            distance;
  };
  std::vector<DiagonalDistance> diagonalDistances;

  for (uint32_t i = 0; i < this->getNcolumns() && i < interactionRadius; i++) {
    for (uint32_t j = i; j < this->getNrows(); j++) {
      auto const dist = qc::NeutralAtomArchitecture::getEuclidianDistance(
          Coordinate(0, 0), Coordinate(i, j));
      if (dist <= interactionRadius) {
        if (dist == 0) {
          continue;
        }
        diagonalDistances.push_back({i, j, dist});
        if (i != j) {
          diagonalDistances.push_back({j, i, dist});
        }
      } else {
        break;
      }
    }
  }
  // sort diagonal distances by distance
  std::sort(diagonalDistances.begin(), diagonalDistances.end(),
            [](DiagonalDistance const& a, DiagonalDistance const& b) {
              return a.distance < b.distance;
            });

  // compute swap distances
  this->swapDistances = SymmetricMatrix(this->getNpositions());

  for (uint32_t coordIndex1 = 0; coordIndex1 < this->getNpositions();
       coordIndex1++) {
    for (uint32_t coordIndex2 = 0; coordIndex2 < coordIndex1; coordIndex2++) {
      auto deltaX = this->getManhattenDistanceX(coordIndex1, coordIndex2);
      auto deltaY = this->getManhattenDistanceY(coordIndex1, coordIndex2);

      // check if one can go diagonal to reduce the swap distance
      uint32_t swapDistance = 0;
      for (auto it = diagonalDistances.rbegin(); it != diagonalDistances.rend();
           ++it) {
        auto const& diagonalDistance = *it;
        while (deltaX >= diagonalDistance.x && deltaY >= diagonalDistance.y) {
          swapDistance += 1;
          deltaX -= diagonalDistance.x;
          deltaY -= diagonalDistance.y;
        }
      }
      // save swap distance in matrix
      this->swapDistances(coordIndex1, coordIndex2) = swapDistance - 1;
      this->swapDistances(coordIndex2, coordIndex1) = swapDistance - 1;
    }
  }
}

void Mapping::swap(Swap swap) {
  auto q1 = swap.first;
  auto q2 = swap.second;
  if (this->isMapped(q1) && this->isMapped(q2)) {
    auto circQ1 = this->getCircQubit(q1);
    auto circQ2 = this->getCircQubit(q2);
    this->removeCircuitQubit(circQ1);
    this->removeCircuitQubit(circQ2);
    this->setCircuitQubit(circQ2, q1);
    this->setCircuitQubit(circQ1, q2);
  } else if (this->isMapped(q1) && !this->isMapped(q2)) {
    this->setCircuitQubit(this->getCircQubit(q1), q2);
    this->removeCircuitQubit(q1);
  } else if (this->isMapped(q2) && !this->isMapped(q1)) {
    this->setCircuitQubit(this->getCircQubit(q2), q1);
    this->removeCircuitQubit(q2);
  } else {
    throw std::runtime_error("Cannot swap unmapped qubits");
  }
}

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
