//
// Created by Ludwig Schmid on 13.10.23.
//

#pragma once

#include "ShuttlingOperation.hpp"
#include "namap/NeutralAtomUtils.hpp"
#include "utils.hpp"

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace qc {
class NeutralAtomArchitecture {
  class Coordinate {
  protected:
    std::uint32_t x;
    std::uint32_t y;

  public:
    Coordinate() = default;
    Coordinate(std::uint32_t x, std::uint32_t y) : x(x), y(y) {}

    [[nodiscard]] inline std::uint32_t getX() const { return x; }
    [[nodiscard]] inline std::uint32_t getY() const { return y; }
    [[nodiscard]] inline std::pair<std::uint32_t, std::uint32_t> getXY() const {
      return std::make_pair(x, y);
    }
    [[nodiscard]] inline fp getEuclidianDistance(const Coordinate& c) const {
      return std::sqrt(std::pow(static_cast<fp>(x) - static_cast<fp>(c.x), 2) +
                       std::pow(static_cast<fp>(y) - static_cast<fp>(c.y), 2));
    }
    [[nodiscard]] inline uint32_t
    getManhattenDistanceX(const Coordinate& c) const {
      if (x > c.x) {
        return x - c.x;
      } else {
        return c.x - x;
      }
    }
    [[nodiscard]] inline uint32_t
    getManhattenDistanceY(const Coordinate& c) const {
      if (y > c.y) {
        return y - c.y;
      } else {
        return c.y - y;
      }
    }
  };

  class Properties {
  protected:
    std::uint16_t nRows;
    std::uint16_t nColumns;
    std::uint16_t nAods;
    std::uint16_t nAodCoordinates;
    fp            interQubitDistance;
    fp            interactionRadius;
    fp            blockingFactor;

  public:
    Properties() = default;
    Properties(std::uint16_t nRows, std::uint16_t nColumns, std::uint16_t nAods,
               std::uint16_t nAodCoordinates, fp interQubitDistance,
               fp interactionRadius, fp blockingFactor)
        : nRows(nRows), nColumns(nColumns), nAods(nAods),
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
    [[nodiscard]] inline fp getInterQubitDistance() const {
      return interQubitDistance;
    }
    [[nodiscard]] inline fp getInteractionRadius() const {
      return interactionRadius;
    }
    [[nodiscard]] inline fp getBlockingFactor() const { return blockingFactor; }
  };

  struct Parameters {
    std::uint32_t               nQubits;
    std::map<OpType, fp>        gateTimes;
    std::map<OpType, fp>        gateAverageFidelities;
    std::map<ShuttlingType, fp> shuttlingTimes;
    std::map<ShuttlingType, fp> shuttlingAverageFidelities;
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
  Properties              properties{};
  Parameters              parameters;
  std::vector<Coordinate> coordinates;
  SymmetricMatrix         swapDistances;

  void createCoordinates();
  void computeSwapDistances(fp interactionRadius);

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
  [[nodiscard]] inline std::uint32_t getNqubits() const {
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
  [[nodiscard]] inline fp getSwapDistance(std::uint32_t idx1,
                                          std::uint32_t idx2) const {
    return swapDistances(idx1, idx2);
  }
  [[nodiscard]] inline fp getSwapDistance(const Coordinate& c1,
                                          const Coordinate& c2) const {
    return swapDistances(c1.getX() + c1.getY() * properties.getNcolumns(),
                         c2.getX() + c2.getY() * properties.getNcolumns());
  }

  [[nodiscard]] inline fp getGateTime(OpType opType) const {
    return parameters.gateTimes.at(opType);
  }
  [[nodiscard]] inline fp getGateAverageFidelity(OpType opType) const {
    return parameters.gateAverageFidelities.at(opType);
  }
  [[nodiscard]] inline fp getShuttlingTime(ShuttlingType shuttlingType) const {
    return parameters.shuttlingTimes.at(shuttlingType);
  }
  [[nodiscard]] inline fp
  getShuttlingAverageFidelity(ShuttlingType shuttlingType) const {
    return parameters.shuttlingAverageFidelities.at(shuttlingType);
  }
  [[nodiscard]] inline fp getDecoherenceTime() const {
    return parameters.decoherenceTimes.getTEff();
  }

  [[nodiscard]] inline std::vector<Coordinate> getCoordinates() const {
    return coordinates;
  }
  [[nodiscard]] inline Coordinate getCoordinate(std::uint32_t idx) const {
    return coordinates[idx];
  }

  [[nodiscard]] inline uint32_t getIndex(const Coordinate& c) {
    return c.getX() + c.getY() * properties.getNcolumns();
  }

  void loadJson(const std::string& filename);

  [[nodiscard]] inline fp getEuclidianDistance(std::uint32_t idx1,
                                               std::uint32_t idx2) const {
    return this->coordinates.at(idx1).getEuclidianDistance(
        this->coordinates.at(idx2));
  }
  [[nodiscard]] inline static fp getEuclidianDistance(const Coordinate& c1,
                                                      const Coordinate& c2) {
    return c1.getEuclidianDistance(c2);
  }
  [[nodiscard]] inline uint32_t
  getManhattenDistanceX(std::uint32_t idx1, std::uint32_t idx2) const {
    return this->coordinates.at(idx1).getManhattenDistanceX(
        this->coordinates.at(idx2));
  }
  [[nodiscard]] inline uint32_t
  getManhattenDistanceY(std::uint32_t idx1, std::uint32_t idx2) const {
    return this->coordinates.at(idx1).getManhattenDistanceY(
        this->coordinates.at(idx2));
  }
  [[nodiscard]] inline uint32_t getManhattenDistanceX(const Coordinate& c1) {
    return c1.getManhattenDistanceX(c1);
  }
  [[nodiscard]] inline uint32_t getManhattenDistanceY(const Coordinate& c1) {
    return c1.getManhattenDistanceY(c1);
  }
};

} // namespace qc
