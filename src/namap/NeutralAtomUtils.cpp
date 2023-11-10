//
// Created by Ludwig Schmid on 16.10.23.
//

#include "namap/NeutralAtomUtils.hpp"

namespace qc {

bool MoveVector::overlap(const qc::MoveVector& other) const {
  // do not consider direction for overlap
  const auto firstStartX  = std::min(xStart, xEnd);
  const auto firstEndX    = std::max(xStart, xEnd);
  const auto secondStartX = std::min(other.xStart, other.xEnd);
  const auto secondEndX   = std::max(other.xStart, other.xEnd);
  const auto firstStartY  = std::min(yStart, yEnd);
  const auto firstEndY    = std::max(yStart, yEnd);
  const auto secondStartY = std::min(other.yStart, other.yEnd);
  const auto secondEndY   = std::max(other.yStart, other.yEnd);

  // need to compute all combinations, as sometimes the start and end x/y points
  // are the same
  auto overlapXFirstStart =
      firstStartX >= secondStartX && firstStartX <= secondEndX;
  auto overlapXFirstEnd = firstEndX >= secondStartX && firstEndX <= secondEndX;
  auto overlapXSecondStart =
      secondStartX >= firstStartX && secondStartX <= firstEndX;
  auto overlapXSecondEnd = secondEndX >= firstStartX && secondEndX <= firstEndX;
  auto overlapYFirstStart =
      firstStartY >= secondStartY && firstStartY <= secondEndY;
  auto overlapYFirstEnd = firstEndY >= secondStartY && firstEndY <= secondEndY;
  auto overlapYSecondStart =
      secondStartY >= firstStartY && secondStartY <= firstEndY;
  auto overlapYSecondEnd = secondEndY >= firstStartY && secondEndY <= firstEndY;

  return (overlapXFirstStart || overlapXFirstEnd || overlapXSecondStart ||
          overlapXSecondEnd || overlapYFirstStart || overlapYFirstEnd ||
          overlapYSecondStart || overlapYSecondEnd);
}

bool MoveVector::include(const qc::MoveVector& other) const {
  const auto firstStartX  = std::min(xStart, xEnd);
  const auto firstEndX    = std::max(xStart, xEnd);
  const auto secondStartX = std::min(other.xStart, other.xEnd);
  const auto secondEndX   = std::max(other.xStart, other.xEnd);
  const auto firstStartY  = std::min(yStart, yEnd);
  const auto firstEndY    = std::max(yStart, yEnd);
  const auto secondStartY = std::min(other.yStart, other.yEnd);
  const auto secondEndY   = std::max(other.yStart, other.yEnd);

  const auto includeX =
      (secondStartX < firstStartX) && (firstEndX < secondEndX);
  const auto includeY =
      (secondStartY < firstStartY) && (firstEndY < secondEndY);

  return includeX || includeY;
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
