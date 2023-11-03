//
// Created by Ludwig Schmid on 16.10.23.
//

#include "namap/NeutralAtomUtils.hpp"

namespace qc {

bool MoveVector::overlap(const qc::MoveVector& other) const {
  auto overlapXStart = xStart <= other.xStart && other.xStart <= xEnd;
  auto overlapXEnd   = xStart <= other.xEnd && other.xEnd <= xEnd;
  auto overlapYStart = yStart <= other.yStart && other.yStart <= yEnd;
  auto overlapYEnd   = yStart <= other.yEnd && other.yEnd <= yEnd;
  return (overlapXStart || overlapXEnd) && (overlapYStart || overlapYEnd);
}

bool MoveVector::include(const qc::MoveVector& other) const {
  auto includeX =
      std::signbit(xStart - other.xStart) != std::signbit(xEnd - other.xEnd);
  auto includeY =
      std::signbit(yStart - other.yStart) != std::signbit(yEnd - other.yEnd);
  return includeX && includeY;
}

void MoveCombs::addMoveComb(const qc::MoveComb& otherMove) {
  for (auto& comb : moveCombs) {
    if (comb == otherMove) {
      comb.cost = std::numeric_limits<fp>::quiet_NaN();
      return;
    }
  }
  moveCombs.push_back(otherMove);
}

void MoveCombs::addMoveCombs(const qc::MoveCombs& otherMoveCombs) {
  for (const auto& otherMove : otherMoveCombs.moveCombs) {
    addMoveComb(otherMove);
  }
}

void MoveCombs::removeAllWithSameStart(const qc::MoveComb& moveComb) {
  moveCombs.erase(std::remove_if(moveCombs.begin(), moveCombs.end(),
                                 [&moveComb](const MoveComb& comb) {
                                   return comb.getFirstMove().first ==
                                          moveComb.getFirstMove().first;
                                 }),
                  moveCombs.end());
}

void MoveCombs::removeAllWithSameEnd(const qc::MoveComb& moveComb) {
  moveCombs.erase(std::remove_if(moveCombs.begin(), moveCombs.end(),
                                 [&moveComb](const MoveComb& comb) {
                                   return comb.getLastMove().second ==
                                          moveComb.getLastMove().second;
                                 }),
                  moveCombs.end());
}

} // namespace qc
