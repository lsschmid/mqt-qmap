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
  this->swapDistances = SwapDistances(
      this->getNpositions(), std::vector<uint32_t>(this->getNpositions(), 0));

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
      this->swapDistances[coordIndex1][coordIndex2] = swapDistance - 1;
      this->swapDistances[coordIndex2][coordIndex1] = swapDistance - 1;
    }
  }
}

} // namespace qc
