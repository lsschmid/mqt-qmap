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

} // namespace qc
