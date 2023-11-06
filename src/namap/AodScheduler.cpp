//
// Created by Ludwig Schmid on 30.10.23.
//

#include "namap/AodScheduler.hpp"

#include "operations/AodOperation.hpp"

namespace qc {

QuantumComputation AodScheduler::schedule(QuantumComputation& qc) {
  initMoveGroups(qc);
  if (moveGroups.empty()) {
    return qc;
  }
  processMoveGroups();

  auto     groupIt = moveGroups.begin();
  uint32_t idx     = 0;
  for (const auto& op : qc) {
    if (groupIt != moveGroups.end() && idx == groupIt->getFirstIdx()) {
      // add move group
      for (auto& aodOp : groupIt->processedOps) {
        qcScheduled.emplace_back(aodOp->clone());
      }
      groupIt++;
    } else {
      qcScheduled.emplace_back(op->clone());
    }
    idx++;
  }

  return qcScheduled;
}

void AodScheduler::initMoveGroups(QuantumComputation& qc) {
  MoveGroup       currentMoveGroup{arch};
  MoveGroup const lastMoveGroup{arch};
  uint32_t        idx = 0;
  for (auto& op : qc) {
    if (op->getType() == OpType::Move) {
      AtomMove const move{op->getTargets()[0], op->getTargets()[1]};
      if (currentMoveGroup.canAdd(move)) {
        currentMoveGroup.add(move, idx);
      } else {
        moveGroups.push_back(std::move(currentMoveGroup));
        currentMoveGroup = MoveGroup{arch};
        currentMoveGroup.add(move, idx);
      }
    }
    idx++;
  }
  if (!currentMoveGroup.moves.empty()) {
    moveGroups.push_back(std::move(currentMoveGroup));
  }
}

bool AodScheduler::MoveGroup::canAdd(const AtomMove& move) {
  // checks if the op can be executed in parallel
  auto moveVector = arch.getVector(move.first, move.second);
  return std::all_of(
      moves.begin(), moves.end(),
      [&moveVector, this](const std::pair<AtomMove, uint32_t> opPair) {
        auto moveGroup = opPair.first;
        auto opVector  = arch.getVector(moveGroup.first, moveGroup.second);
        return parallelCheck(moveVector, opVector);
      });
}

bool AodScheduler::MoveGroup::parallelCheck(const MoveVector& v1,
                                            const MoveVector& v2) {
  if (!v1.overlap(v2)) {
    return true;
  }
  // overlap -> check if same direction
  if (v1.direction != v2.direction) {
    return false;
  }
  // same direction -> check if include
  if (v1.include(v2) || v2.include(v1)) {
    return false;
  }
  return true;
}

void AodScheduler::MoveGroup::add(const AtomMove& move, const uint32_t idx) {
  moves.emplace_back(move, idx);
}

bool AodScheduler::AodActivationHelper::addActivation(const Coordinate& origin,
                                                      const AtomMove&   move,
                                                      MoveVector        v) {
  const auto          x         = origin.getX();
  const auto          y         = origin.getY();
  const auto          signX     = v.direction.getSignX();
  const auto          signY     = v.direction.getSignY();
  const auto          deltaX    = v.xEnd - v.xStart;
  const auto          deltaY    = v.yEnd - v.yStart;
  auto                aodMovesX = getAodMoves(Dimension::X, x);
  auto                aodMovesY = getAodMoves(Dimension::Y, y);
  AodActivation const newActivationPair{
      {x, deltaX, signX}, {y, deltaY, signY}, move};
  if (aodMovesX.empty() && aodMovesY.empty()) {
    // create new activation
    allActivations.push_back(newActivationPair);
    return true;
  }
  if (aodMovesY.empty()) {
    // check if it can be combined with existing activations
    for (auto* aodMove : aodMovesX) {
      if (aodMove->init == x && aodMove->delta == deltaX &&
          aodMove->offset == signX) {
        // combine activations
        mergeActivation(Dimension::X, newActivationPair);
        return true;
      }
    }

    // check if increase in x direction possible
    if (!checkIntermediateSpace(Dimension::X, x, signX)) {
      return false;
    }
    allActivations.push_back(newActivationPair);
    // sort aodMovesX by delta and assign offset
    reAssignOffsets(aodMovesX, signX);
    return true;
  }
  if (aodMovesX.empty()) {
    // check if it can be combined with existing activations
    for (auto* aodMove : aodMovesY) {
      if (aodMove->init == y && aodMove->delta == deltaY &&
          aodMove->offset == signY) {
        // combine activations
        mergeActivation(Dimension::Y, newActivationPair);
        return true;
      }
    }

    // check if increase in y direction possible
    if (!checkIntermediateSpace(Dimension::Y, y, signY)) {
      return false;
    }
    allActivations.push_back(newActivationPair);
    // sort aodMovesY by delta and assign offset
    reAssignOffsets(aodMovesY, signY);
    return true;
  }
  // check if it can be combined with existing activations
  bool combinedX = false;
  bool combinedY = false;

  for (auto* aodMove : aodMovesX) {
    if (aodMove->init == x && aodMove->delta == deltaX &&
        aodMove->offset == signX) {
      // combine activations
      mergeActivation(Dimension::X, newActivationPair);
      combinedX = true;
    }
  }
  for (auto* aodMove : aodMovesY) {
    if (aodMove->init == y && aodMove->delta == deltaY &&
        aodMove->offset == signY) {
      // combine activations
      mergeActivation(Dimension::Y, newActivationPair);
      combinedY = true;
    }
  }
  if (combinedX && combinedY) {
    return true;
  }
  if (!checkIntermediateSpace(Dimension::X, x, signX) ||
      !checkIntermediateSpace(Dimension::Y, y, signY)) {
    return false;
  }
  allActivations.push_back(newActivationPair);
  reAssignOffsets(aodMovesX, signX);
  reAssignOffsets(aodMovesY, signY);
  return true;
}

void AodScheduler::AodActivationHelper::reAssignOffsets(
    std::vector<AodMove*>& aodMoves, int32_t sign) {
  std::sort(aodMoves.begin(), aodMoves.end(),
            [](const AodMove* a, const AodMove* b) {
              return std::abs(a->delta) < std::abs(b->delta);
            });
  int32_t offset = sign;
  for (auto* aodMove : aodMoves) {
    if (std::signbit(aodMove->delta) == std::signbit(sign)) {
      aodMove->offset = offset;
    }
    offset += sign;
  }
}

void AodScheduler::processMoveGroups() {
  // convert the moves from MoveGroup to AodOperations
  for (auto groupIt = moveGroups.begin(); groupIt != moveGroups.end();
       ++groupIt) {
    auto&               moveGroup = *groupIt;
    AodActivationHelper aodActivationHelper{arch};
    AodActivationHelper aodDeactivationHelper{arch};
    MoveGroup           possibleNewMoveGroup{arch};
    for (auto& movePair : moveGroup.moves) {
      auto& move   = movePair.first;
      auto  idx    = movePair.second;
      auto  origin = arch.getCoordinate(move.first);
      auto  v      = arch.getVector(move.first, move.second);
      if (!aodActivationHelper.addActivation(origin, move, v)) {
        // move could not be added as not sufficient intermediate levels
        // add new move group and add move to it
        possibleNewMoveGroup.add(move, idx);
      }
    }
    if (!possibleNewMoveGroup.moves.empty()) {
      groupIt = moveGroups.insert(groupIt + 1, std::move(possibleNewMoveGroup));
    }
    groupIt->processedOps = aodActivationHelper.getAodOperations();
  }
}

std::vector<AodScheduler::AodActivationHelper::AodMove*>
AodScheduler::AodActivationHelper::getAodMoves(Dimension dim,
                                               uint32_t  x) const {
  std::vector<AodMove*> aodMoves;
  for (const auto& activation : allActivations) {
    for (auto* aodMove : activation.getActivates(dim)) {
      if (aodMove->init == x) {
        aodMoves.push_back(aodMove);
      }
    }
  }
  return aodMoves;
}

uint32_t AodScheduler::AodActivationHelper::getMaxOffset(Dimension dim,
                                                         uint32_t  x,
                                                         int32_t   sign) const {
  auto aodMoves = getAodMoves(dim, x);
  if (aodMoves.empty()) {
    return 0;
  }
  uint32_t maxOffset = 0;
  for (const auto& aodMove : aodMoves) {
    auto offset = aodMove->offset;
    if (offset >= sign) {
      maxOffset = std::max(maxOffset, static_cast<uint32_t>(offset));
    } else {
      maxOffset = std::max(maxOffset, static_cast<uint32_t>(-offset));
    }
  }
  return maxOffset;
}

bool AodScheduler::AodActivationHelper::checkIntermediateSpace(
    Dimension dim, uint32_t x, int32_t sign) const {
  uint32_t neighborX = x;
  if (sign > 0) {
    neighborX += 1;
  } else {
    neighborX -= 1;
  }
  auto aodMoves         = getAodMoves(dim, x);
  auto aodMovesNeighbor = getAodMoves(dim, neighborX);
  if (aodMoves.empty() && aodMovesNeighbor.empty()) {
    return true;
  }
  if (aodMoves.empty()) {
    return getMaxOffset(dim, neighborX, sign) < nIntermediateLevels;
  }
  if (aodMovesNeighbor.empty()) {
    return getMaxOffset(dim, x, sign) < nIntermediateLevels;
  }
  return getMaxOffset(dim, x, sign) + getMaxOffset(dim, neighborX, sign) <
         nIntermediateLevels;
}

void AodScheduler::AodActivationHelper::mergeActivation(
    Dimension dim, const AodActivationHelper::AodActivation& activation) {
  // merge activations
  for (auto& activationCurrent : allActivations) {
    auto activates = activation.getActivates(dim);
    for (auto* aodMove : activates) {
      if (aodMove->init == activation.getActivates(dim)[0]->init &&
          aodMove->delta == activation.getActivates(dim)[0]->delta &&
          aodMove->offset == activation.getActivates(dim)[0]->offset) {
        activationCurrent.moves.push_back(activation.moves[0]);
        Dimension const otherDim =
            dim == Dimension::X ? Dimension::Y : Dimension::X;
        auto* addedAodMove = activation.getActivates(otherDim)[0];
        activationCurrent.addAodMove(otherDim, *addedAodMove);
        // sort aodMovesY by delta and assign offset
        auto aodMoves = getAodMoves(otherDim, addedAodMove->init);
        reAssignOffsets(aodMoves, addedAodMove->offset);
        return;
      }
    }
  }
}

std::pair<OpPointer, OpPointer>
AodScheduler::AodActivationHelper::getAodOperation(
    const AodActivationHelper::AodActivation& activation,
    const NeutralAtomArchitecture&            arch) {
  std::vector<CoordIndex> qubitsActivation;
  qubitsActivation.reserve(activation.moves.size());
  for (const auto& move : activation.moves) {
    qubitsActivation.push_back(move.first);
  }
  std::vector<CoordIndex> qubitsMove;
  qubitsMove.reserve(activation.moves.size() * 2);
  for (const auto& move : activation.moves) {
    qubitsMove.push_back(move.first);
    qubitsMove.push_back(move.second);
  }

  std::vector<std::tuple<uint32_t, fp, fp>> initOperations;
  std::vector<std::tuple<uint32_t, fp, fp>> offsetOperations;

  auto d      = arch.getInterQubitDistance();
  auto interD = arch.getInterQubitDistance() / arch.getNAodIntermediateLevels();

  for (const auto& aodMove : activation.activateXs) {
    initOperations.emplace_back(0, static_cast<fp>(aodMove->init) * d,
                                static_cast<fp>(aodMove->init) * d);
    offsetOperations.emplace_back(0, static_cast<fp>(aodMove->init) * d,
                                  static_cast<fp>(aodMove->init) * d +
                                      static_cast<fp>(aodMove->offset) *
                                          interD);
  }
  for (const auto& aodMove : activation.activateYs) {
    initOperations.emplace_back(1, static_cast<fp>(aodMove->init) * d,
                                static_cast<fp>(aodMove->init) * d);
    offsetOperations.emplace_back(1, static_cast<fp>(aodMove->init) * d,
                                  static_cast<fp>(aodMove->init) * d +
                                      static_cast<fp>(aodMove->offset) *
                                          interD);
  }

  auto initOp = std::make_unique<AodOperation>(
      OpType::AodActivate, qubitsActivation, initOperations);
  auto offsetOp = std::make_unique<AodOperation>(OpType::AodMove, qubitsMove,
                                                 offsetOperations);
  return std::make_pair(std::move(initOp), std::move(offsetOp));
}

std::vector<OpPointer>
AodScheduler::AodActivationHelper::getAodOperations() const {
  std::vector<OpPointer> aodOperations;
  for (const auto& activation : allActivations) {
    auto [initOp, offsetOp] = getAodOperation(activation, arch);
    aodOperations.emplace_back(std::move(initOp));
    aodOperations.emplace_back(std::move(offsetOp));
  }
  return aodOperations;
}

} // namespace qc
