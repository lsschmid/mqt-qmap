//
// Created by Ludwig Schmid on 16.10.23.
//

#pragma once

#include "Definitions.hpp"
#include "utils.hpp"

namespace qc {

// Symmetric matrix class with same number of rows and columns that allows
// access by row and column but uses less memory than a full matrix

class SymmetricMatrix {
private:
  std::vector<std::vector<fp>> data;
  uint32_t                     size = 0;

public:
  SymmetricMatrix() = default;
  explicit SymmetricMatrix(uint32_t size) : size(size) {
    data.resize(size);
    for (uint32_t i = 0; i < size; ++i) {
      data[i].resize(i + 1);
    }
  }

  SymmetricMatrix(uint32_t size, fp value) : size(size) {
    data.resize(size);
    for (uint32_t i = 0; i < size; ++i) {
      data[i].resize(i + 1, value);
    }
  }

  inline fp& operator()(uint32_t row, uint32_t col) {
    if (row < col) {
      return data[col][row];
    }
    return data[row][col];
  }

  [[nodiscard]] inline fp operator()(uint32_t row, uint32_t col) const {
    if (row < col) {
      return data[col][row];
    }
    return data[row][col];
  }

  [[nodiscard]] inline uint32_t getSize() const { return size; }
};

enum InitialCoordinateMapping { Trivial, Random };
enum InitialMapping { Identity };

struct MoveVector {
  struct Direction {
    bool x;
    bool y;

    Direction(bool x, bool y) : x(x), y(y) {}
    Direction(fp deltaX, fp deltaY) : x(deltaX >= 0), y(deltaY >= 0) {}

    [[nodiscard]] inline bool operator==(const Direction& other) const {
      return x == other.x && y == other.y;
    }
    [[nodiscard]] inline bool operator!=(const Direction& other) const {
      return !(*this == other);
    }
  };

  fp        xStart;
  fp        yStart;
  fp        xEnd;
  fp        yEnd;
  Direction direction;

  MoveVector(fp xStart, fp yStart, fp xEnd, fp yEnd)
      : xStart(xStart), yStart(yStart), xEnd(xEnd), yEnd(yEnd),
        direction(xEnd - xStart, yEnd - yStart) {}

  [[nodiscard]] inline bool sameDirection(const MoveVector& other) const {
    return direction == other.direction;
  }
  [[nodiscard]] inline fp getLength() const {
    return std::sqrt(std::pow(xEnd - xStart, 2) + std::pow(yEnd - yStart, 2));
  }
  [[nodiscard]] bool overlap(const MoveVector& other) const;
  [[nodiscard]] bool include(const MoveVector& other) const;
};

} // namespace qc
