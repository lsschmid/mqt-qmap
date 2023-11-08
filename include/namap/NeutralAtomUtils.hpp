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
  [[nodiscard]] inline int32_t getSignX() const { return x ? 1 : -1; }
  [[nodiscard]] inline int32_t getSignY() const { return y ? 1 : -1; }
};

struct MoveVector {
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

  MoveComb(std::vector<AtomMove> moves, fp cost)
      : moves(std::move(moves)), cost(cost) {}
  MoveComb(AtomMove move, fp cost)
      : moves(std::vector<AtomMove>{std::move(move)}), cost(cost) {}

  MoveComb() = default;
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

  inline void append(AtomMove& addMove) {
    moves.push_back(addMove);
    cost = std::numeric_limits<fp>::quiet_NaN();
  }
  inline void append(const MoveComb& addMoveComb) {
    moves.insert(moves.end(), addMoveComb.moves.begin(),
                 addMoveComb.moves.end());
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

class Coordinate {
protected:
  std::uint32_t x;
  std::uint32_t y;

public:
  Coordinate() = default;
  Coordinate(CoordIndex x, CoordIndex y) : x(x), y(y) {}

  [[nodiscard]] inline CoordIndex getX() const { return x; }
  [[nodiscard]] inline CoordIndex getY() const { return y; }
  [[nodiscard]] inline fp         getXfp() const { return static_cast<fp>(x); }
  [[nodiscard]] inline fp         getYfp() const { return static_cast<fp>(y); }
  [[nodiscard]] inline std::pair<CoordIndex, CoordIndex> getXY() const {
    return std::make_pair(x, y);
  }
  [[nodiscard]] inline fp getEuclidianDistance(const Coordinate& c) const {
    return std::sqrt(std::pow(static_cast<fp>(x) - static_cast<fp>(c.x), 2) +
                     std::pow(static_cast<fp>(y) - static_cast<fp>(c.y), 2));
  }
  [[nodiscard]] static inline bool sameX(const Coordinate& c1,
                                         const Coordinate& c2) {
    return c1.x == c2.x;
  }
  [[nodiscard]] static inline bool sameY(const Coordinate& c1,
                                         const Coordinate& c2) {
    return c1.y == c2.y;
  }
  [[nodiscard]] static inline bool sameXorY(const Coordinate& c1,
                                            const Coordinate& c2) {
    return sameX(c1, c2) || sameY(c1, c2);
  }

  [[nodiscard]] inline CoordIndex
  getManhattenDistanceX(const Coordinate& c) const {
    if (x > c.x) {
      return x - c.x;
    } else {
      return c.x - x;
    }
  }
  [[nodiscard]] inline CoordIndex
  getManhattenDistanceY(const Coordinate& c) const {
    if (y > c.y) {
      return y - c.y;
    } else {
      return c.y - y;
    }
  }
};

} // namespace qc
