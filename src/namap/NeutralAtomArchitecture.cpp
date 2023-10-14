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
        jsonDataProperties["interQubitDistance"]);

    // Load parameters
    const nlohmann::json jsonDataParameters = jsonData["parameters"];
    this->parameters                        = Parameters();
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

[[nodiscard]] fp
NeutralAtomArchitecture::getDistance(std::uint32_t idx1,
                                     std::uint32_t idx2) const {
  return this->coordinates.at(idx1).getDistance(this->coordinates.at(idx2));
}

[[nodiscard]] fp NeutralAtomArchitecture::getDistance(
    const qc::NeutralAtomArchitecture::Coordinate& c1,
    const qc::NeutralAtomArchitecture::Coordinate& c2) {
  return c1.getDistance(c2);
}

} // namespace qc
