//
// Created by Ludwig Schmid on 13.10.23.
//

#pragma once

#include "ShuttlingOperation.hpp"

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
    [[nodiscard]] inline fp getDistance(const Coordinate& c) const {
      return std::sqrt(std::pow(static_cast<fp>(x) - static_cast<fp>(c.x), 2) +
                       std::pow(static_cast<fp>(y) - static_cast<fp>(c.y), 2));
    }
  };

  class Properties {
  protected:
    std::uint16_t nRows;
    std::uint16_t nColumns;
    std::uint16_t nAods;
    std::uint16_t nAodCoordinates;
    fp            interQubitDistance;

  public:
    Properties() = default;
    Properties(std::uint16_t nRows, std::uint16_t nColumns, std::uint16_t nAods,
               std::uint16_t nAodCoordinates, fp interQubitDistance)
        : nRows(nRows), nColumns(nColumns), nAods(nAods),
          nAodCoordinates(nAodCoordinates),
          interQubitDistance(interQubitDistance) {}
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
  };

  struct Parameters {
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

  void createCoordinates();

public:
  std::string name;

  NeutralAtomArchitecture() = default;
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
  [[nodiscard]] inline fp getInterQubitDistance() const {
    return properties.getInterQubitDistance();
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

  void loadJson(const std::string& filename);

  [[nodiscard]] fp getDistance(std::uint32_t idx1, std::uint32_t idx2) const;
  [[nodiscard]] static fp getDistance(const Coordinate& c1,
                                      const Coordinate& c2);
};
} // namespace qc
