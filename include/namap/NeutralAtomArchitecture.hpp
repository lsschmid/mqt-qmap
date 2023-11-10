//
// Created by Ludwig Schmid on 13.10.23.
//

#pragma once

#include "namap/NeutralAtomDefinitions.hpp"
#include "namap/NeutralAtomUtils.hpp"
#include "utils.hpp"

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace qc {
class NeutralAtomArchitecture {
  class Properties {
  protected:
    std::uint16_t nRows;
    std::uint16_t nColumns;
    std::uint16_t nAods;
    std::uint16_t nAodIntermediateLevels;
    std::uint16_t nAodCoordinates;
    fp            interQubitDistance;
    fp            interactionRadius;
    fp            blockingFactor;

  public:
    Properties() = default;
    Properties(std::uint16_t nRows, std::uint16_t nColumns, std::uint16_t nAods,
               std::uint16_t nAodCoordinates, fp interQubitDistance,
               fp interactionRadius, fp blockingFactor, fp aodDistance)
        : nRows(nRows), nColumns(nColumns), nAods(nAods),
          nAodIntermediateLevels(
              static_cast<uint16_t>(interQubitDistance / aodDistance)),
          nAodCoordinates(nAodCoordinates),
          interQubitDistance(interQubitDistance),
          interactionRadius(interactionRadius), blockingFactor(blockingFactor) {
    }
    [[nodiscard]] inline std::uint16_t getNpositions() const {
      return nRows * nColumns;
    }
    [[nodiscard]] inline std::uint16_t getNrows() const { return nRows; }
    [[nodiscard]] inline std::uint16_t getNcolumns() const { return nColumns; }
    [[nodiscard]] inline std::uint16_t getNAods() const { return nAods; }
    [[nodiscard]] inline std::uint16_t getNAodCoordinates() const {
      return nAodCoordinates;
    }
    [[nodiscard]] inline std::uint16_t getNAodIntermediateLevels() const {
      return nAodIntermediateLevels;
    }
    [[nodiscard]] inline fp getInterQubitDistance() const {
      return interQubitDistance;
    }
    [[nodiscard]] inline fp getInteractionRadius() const {
      return interactionRadius;
    }
    [[nodiscard]] inline fp getBlockingFactor() const { return blockingFactor; }
  };

  struct Parameters {
    CoordIndex                nQubits;
    std::map<std::string, fp> gateTimes;
    std::map<std::string, fp> gateAverageFidelities;
    std::map<OpType, fp>      shuttlingTimes;
    std::map<OpType, fp>      shuttlingAverageFidelities;
    class DecoherenceTimes {
    protected:
      fp tEff;

    public:
      fp t1;
      fp t2;
      DecoherenceTimes() = default;
      inline DecoherenceTimes(fp t1, fp t2)
          : tEff(t1 * t2 / (t1 + t2)), t1(t1), t2(t2) {}
      [[nodiscard]] inline fp getTEff() const { return tEff; }
    };
    DecoherenceTimes decoherenceTimes;
  };

protected:
  Properties                        properties{};
  Parameters                        parameters;
  std::vector<Coordinate>           coordinates;
  SymmetricMatrix                   swapDistances;
  std::vector<std::set<CoordIndex>> nearbyCoordinates;

  void createCoordinates();
  void computeSwapDistances(fp interactionRadius);
  void computeNearbyCoordinates();

public:
  std::string name;

  NeutralAtomArchitecture() = delete;
  // create from JSON file
  explicit NeutralAtomArchitecture(const std::string& filename);

  [[nodiscard]] inline std::uint16_t getNrows() const {
    return properties.getNrows();
  }
  [[nodiscard]] inline std::uint16_t getNcolumns() const {
    return properties.getNcolumns();
  }
  [[nodiscard]] inline std::uint16_t getNpositions() const {
    return properties.getNpositions();
  }
  [[nodiscard]] inline std::uint16_t getNAods() const {
    return properties.getNAods();
  }
  [[nodiscard]] inline std::uint16_t getNAodCoordinates() const {
    return properties.getNAodCoordinates();
  }
  [[nodiscard]] inline CoordIndex getNqubits() const {
    return parameters.nQubits;
  }
  [[nodiscard]] inline fp getInterQubitDistance() const {
    return properties.getInterQubitDistance();
  }

  [[nodiscard]] inline fp getInteractionRadius() const {
    return properties.getInteractionRadius();
  }
  [[nodiscard]] inline fp getBlockingFactor() const {
    return properties.getBlockingFactor();
  }
  [[nodiscard]] inline fp getSwapDistance(CoordIndex idx1,
                                          CoordIndex idx2) const {
    return swapDistances(idx1, idx2);
  }
  [[nodiscard]] inline fp getSwapDistance(const Coordinate& c1,
                                          const Coordinate& c2) const {
    return swapDistances(c1.getX() + c1.getY() * properties.getNcolumns(),
                         c2.getX() + c2.getY() * properties.getNcolumns());
  }

  [[nodiscard]] inline uint16_t getNAodIntermediateLevels() const {
    return properties.getNAodIntermediateLevels();
  }

  [[nodiscard]] fp getOpTime(const Operation* op) const;
  [[nodiscard]] fp getOpFidelity(const Operation* op) const;
  [[nodiscard]] std::set<CoordIndex>
  getBlockedCoordIndices(const Operation* op) const;

  [[nodiscard]] inline fp getGateTime(std::string s) const {
    if (parameters.gateTimes.find(s) == parameters.gateTimes.end()) {
      std::cout << "Gate time for " << s << " not found\n"
                << "Returning default value\n";
      return parameters.gateTimes.at("none");
    }
    return parameters.gateTimes.at(s);
  }
  [[nodiscard]] inline fp getGateAverageFidelity(std::string s) const {
    if (parameters.gateAverageFidelities.find(s) ==
        parameters.gateAverageFidelities.end()) {
      std::cout << "Gate average fidelity for " << s << " not found\n"
                << "Returning default value\n";
      return parameters.gateAverageFidelities.at("none");
    }
    return parameters.gateAverageFidelities.at(s);
  }
  [[nodiscard]] inline fp getShuttlingTime(OpType shuttlingType) const {
    return parameters.shuttlingTimes.at(shuttlingType);
  }
  [[nodiscard]] inline fp
  getShuttlingAverageFidelity(OpType shuttlingType) const {
    return parameters.shuttlingAverageFidelities.at(shuttlingType);
  }
  [[nodiscard]] inline fp getDecoherenceTime() const {
    return parameters.decoherenceTimes.getTEff();
  }

  [[nodiscard]] inline std::vector<Coordinate> getCoordinates() const {
    return coordinates;
  }
  [[nodiscard]] inline Coordinate getCoordinate(CoordIndex idx) const {
    return coordinates[idx];
  }

  [[nodiscard]] inline CoordIndex getIndex(const Coordinate& c) {
    return c.getX() + c.getY() * properties.getNcolumns();
  }

  void loadJson(const std::string& filename);

  [[nodiscard]] inline fp getEuclidianDistance(CoordIndex idx1,
                                               CoordIndex idx2) const {
    return this->coordinates.at(idx1).getEuclidianDistance(
        this->coordinates.at(idx2));
  }
  [[nodiscard]] inline static fp getEuclidianDistance(const Coordinate& c1,
                                                      const Coordinate& c2) {
    return c1.getEuclidianDistance(c2);
  }
  [[nodiscard]] inline CoordIndex getManhattenDistanceX(CoordIndex idx1,
                                                        CoordIndex idx2) const {
    return this->coordinates.at(idx1).getManhattenDistanceX(
        this->coordinates.at(idx2));
  }
  [[nodiscard]] inline CoordIndex getManhattenDistanceY(CoordIndex idx1,
                                                        CoordIndex idx2) const {
    return this->coordinates.at(idx1).getManhattenDistanceY(
        this->coordinates.at(idx2));
  }
  [[nodiscard]] inline CoordIndex getManhattenDistanceX(const Coordinate& c1) {
    return c1.getManhattenDistanceX(c1);
  }
  [[nodiscard]] inline CoordIndex getManhattenDistanceY(const Coordinate& c1) {
    return c1.getManhattenDistanceY(c1);
  }

  [[nodiscard]] inline std::set<CoordIndex>
  getNearbyCoordinates(CoordIndex idx) const {
    return nearbyCoordinates[idx];
  }
  [[nodiscard]] inline MoveVector getVector(CoordIndex idx1, CoordIndex idx2) {
    return {this->coordinates[idx1].getXfp(), this->coordinates[idx1].getYfp(),
            this->coordinates[idx2].getXfp(), this->coordinates[idx2].getYfp()};
  }

  [[nodiscard]] std::vector<CoordIndex> getNN(CoordIndex idx) const;

  [[nodiscard]] inline fp getVectorShuttlingTime(const MoveVector& v) const {
    return v.getLength() * this->getInterQubitDistance() /
           this->getShuttlingTime(OpType::Move);
  }
};

} // namespace qc
