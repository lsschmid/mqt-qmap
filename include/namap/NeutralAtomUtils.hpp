//
// Created by Ludwig Schmid on 16.10.23.
//

#pragma once

#include "Definitions.hpp"
#include "namap/NeutralAtomDefinitions.hpp"
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

struct MoveComb {
  std::vector<AtomMove> moves;
  fp                    cost = std::numeric_limits<fp>::quiet_NaN();

  MoveComb(std::vector<AtomMove> moves, fp weight, fp cost)
      : moves(std::move(moves)), cost(cost) {}
  MoveComb(AtomMove, fp weight, fp cost)
      : moves(std::vector<AtomMove>{std::move(moves)}), cost(cost) {}

  MoveComb(std::vector<AtomMove> moves) : moves(std::move(moves)) {}
  MoveComb(AtomMove move) : moves(std::vector<AtomMove>{std::move(move)}) {}

  [[nodiscard]] AtomMove inline getFirstMove() const { return *moves.begin(); }
  [[nodiscard]] AtomMove inline getLastMove() const { return *moves.rbegin(); }
  // implement == operator for AtomMove
  [[nodiscard]] inline bool operator==(const MoveComb& other) const {
    return moves == other.moves;
  }
  [[nodiscard]] inline bool operator!=(const MoveComb& other) const {
    return !(*this == other);
  }

  [[nodiscard]] inline void append(AtomMove& addMove) {
    moves.push_back(addMove);
    cost = std::numeric_limits<fp>::quiet_NaN();
  }
  [[nodiscard]] inline size_t size() const { return moves.size(); }
  [[nodiscard]] inline bool   empty() const { return moves.empty(); }

  [[nodiscard]] inline bool containsCoord(CoordIndex idx) {
    return std::any_of(moves.begin(), moves.end(), [idx](const AtomMove& move) {
      return move.first == idx || move.second == idx;
    });
  }
};

struct MoveCombs {
  std::vector<MoveComb> moveCombs;

  MoveCombs() = default;
  MoveCombs(std::vector<MoveComb> moveCombs) : moveCombs(moveCombs) {}

  [[nodiscard]] bool inline empty() const { return moveCombs.empty(); }
  [[nodiscard]] size_t inline size() const { return moveCombs.size(); }

  // define iterators that iterate over the moveCombs vector
  typedef std::vector<MoveComb>::iterator       iterator;
  typedef std::vector<MoveComb>::const_iterator const_iterator;
  iterator       begin() { return moveCombs.begin(); }
  iterator       end() { return moveCombs.end(); }
  const_iterator begin() const { return moveCombs.begin(); }
  const_iterator end() const { return moveCombs.end(); }

  void addMoveComb(const MoveComb& moveComb);
  void addMoveCombs(const MoveCombs& otherMoveCombs);
  void removeAllWithSameStart(const MoveComb& moveComb);
  void removeAllWithSameEnd(const MoveComb& moveComb);
};

} // namespace qc
