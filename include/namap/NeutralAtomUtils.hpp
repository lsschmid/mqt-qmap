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

enum InitialCoordinateMapping { Trivial };
enum InitialMapping { Identity };

} // namespace qc
