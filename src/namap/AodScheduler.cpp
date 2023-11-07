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
      for (auto& aodOp : groupIt->processedOpsInit) {
        qcScheduled.emplace_back(aodOp->clone());
      }
      groupIt++;
    } else if (op->getType() != OpType::Move) {
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
    } else if (op->getNqubits() > 1) {
      for (const auto& qubit : op->getUsedQubits()) {
        if (std::find(currentMoveGroup.targetQubits.begin(),
                      currentMoveGroup.targetQubits.end(),
                      qubit) == currentMoveGroup.targetQubits.end()) {
          currentMoveGroup.qubitsUsedByGates.push_back(qubit);
        }
      }
    }
    idx++;
  }
  if (!currentMoveGroup.moves.empty()) {
    moveGroups.push_back(std::move(currentMoveGroup));
  }
}

bool AodScheduler::MoveGroup::canAdd(const AtomMove& move) {
  // if move wourld move a qubit that is used by a gate in this move group
  // return false
  if (std::find(qubitsUsedByGates.begin(), qubitsUsedByGates.end(),
                move.first) != qubitsUsedByGates.end()) {
    return false;
  }
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
  targetQubits.push_back(move.second);
}

void AodScheduler::AodActivationHelper::addActivation(
    std::pair<ActivationMerge, ActivationMerge> merge, const Coordinate& origin,
    const AtomMove& move, MoveVector v) {
  const auto x      = origin.getX();
  const auto y      = origin.getY();
  const auto signX  = v.direction.getSignX();
  const auto signY  = v.direction.getSignY();
  const auto deltaX = v.xEnd - v.xStart;
  const auto deltaY = v.yEnd - v.yStart;
  //  AodActivation const newActivationPair{
  //      {x, deltaX, signX}, {y, deltaY, signY}, move};
  //  allActivations.push_back(newActivationPair);
  //  mergeActivation(Dimension::X, newActivationPair);
  //  reAssignOffsets(aodMovesX, signX);
  auto aodMovesX = getAodMoves(Dimension::X, x);
  auto aodMovesY = getAodMoves(Dimension::Y, y);

  switch (merge.first) {
  case ActivationMerge::Trivial:
    switch (merge.second) {
    case ActivationMerge::Trivial:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      break;
    case ActivationMerge::Merge:
      allActivations.emplace_back(
          AodActivation{Dimension::X, {x, deltaX, signX}, move});
      mergeActivation(Dimension::Y,
                      AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      break;
    case ActivationMerge::Append:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      aodMovesY = getAodMoves(Dimension::Y, y);
      reAssignOffsets(aodMovesY, signY);
      break;
    default:
      break;
    }
    break;
  case ActivationMerge::Merge:
    switch (merge.second) {
    case ActivationMerge::Trivial:
      allActivations.emplace_back(
          AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      mergeActivation(Dimension::X,
                      AodActivation{Dimension::X, {x, deltaX, signX}, move});
      break;
    case ActivationMerge::Merge:
      mergeActivation(Dimension::X,
                      AodActivation{Dimension::X, {x, deltaX, signX}, move});
      mergeActivation(Dimension::Y,
                      AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      break;
    case ActivationMerge::Append:
      mergeActivation(Dimension::X,
                      AodActivation{Dimension::X, {x, deltaX, signX}, move});
      allActivations.emplace_back(
          AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      aodMovesY = getAodMoves(Dimension::Y, y);
      reAssignOffsets(aodMovesY, signY);
      break;
    default:
      break;
    }
    break;
  case ActivationMerge::Append:
    switch (merge.second) {
    case ActivationMerge::Trivial:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      aodMovesX = getAodMoves(Dimension::X, x);
      reAssignOffsets(aodMovesX, signX);
      break;
    case ActivationMerge::Merge:
      allActivations.emplace_back(
          AodActivation{Dimension::X, {x, deltaX, signX}, move});
      aodMovesX = getAodMoves(Dimension::X, x);
      reAssignOffsets(aodMovesX, signX);
      mergeActivation(Dimension::Y,
                      AodActivation{Dimension::Y, {y, deltaY, signY}, move});
      break;
    case ActivationMerge::Append:
      allActivations.emplace_back(
          AodActivation{{x, deltaX, signX}, {y, deltaY, signY}, move});
      aodMovesX = getAodMoves(Dimension::X, x);
      reAssignOffsets(aodMovesX, signX);
      aodMovesY = getAodMoves(Dimension::Y, y);
      reAssignOffsets(aodMovesY, signY);
      break;
    default:
      break;
    }
    break;
  default:
    break;
  }
}

std::pair<ActivationMerge, ActivationMerge>
AodScheduler::AodActivationHelper::canAddActivation(
    const qc::Coordinate& origin, const qc::AtomMove& move,
    qc::MoveVector v) const {
  auto aodMovesX = getAodMoves(Dimension::X, origin.getX());
  auto aodMovesY = getAodMoves(Dimension::Y, origin.getY());

  auto const canX = canMergeActivation(Dimension::X, origin, move, v);
  auto const canY = canMergeActivation(Dimension::Y, origin, move, v);
  return std::make_pair(canX, canY);
}

ActivationMerge AodScheduler::AodActivationHelper::canMergeActivation(
    Dimension dim, const Coordinate& origin, const AtomMove& move,
    MoveVector v) const {
  auto x = dim == Dimension::X ? origin.getX() : origin.getY();
  auto sign =
      dim == Dimension::X ? v.direction.getSignX() : v.direction.getSignY();
  auto delta    = dim == Dimension::X ? v.xEnd - v.xStart : v.yEnd - v.yStart;
  auto aodMoves = getAodMoves(dim, x);
  if (aodMoves.empty()) {
    // return false as no merge is required
    return ActivationMerge::Trivial;
  }
  // check if it can be combined with existing activations
  for (auto* aodMove : aodMoves) {
    if (aodMove->init == x && aodMove->delta == delta &&
        aodMove->offset == sign) {
      // combine activations
      return ActivationMerge::Merge;
    }
  }

  // check if increase in x direction possible
  if (checkIntermediateSpace(dim, x, sign)) {
    return ActivationMerge::Append;
  }
  return ActivationMerge::Impossible;
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
      auto& move     = movePair.first;
      auto  idx      = movePair.second;
      auto  origin   = arch.getCoordinate(move.first);
      auto  v        = arch.getVector(move.first, move.second);
      auto  canAddXY = aodActivationHelper.canAddActivation(origin, move, v);
      if (canAddXY.first == ActivationMerge::Impossible ||
          canAddXY.second == ActivationMerge::Impossible) {
        // move could not be added as not sufficient intermediate levels
        // add new move group and add move to it
        possibleNewMoveGroup.add(move, idx);
      } else {
        aodActivationHelper.addActivation(canAddXY, origin, move, v);
      }
    }
    if (!possibleNewMoveGroup.moves.empty()) {
      groupIt = moveGroups.insert(groupIt + 1, std::move(possibleNewMoveGroup));
    }
    groupIt->processedOpsInit = aodActivationHelper.getAodOperations();
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
    return getMaxOffset(dim, neighborX, sign) <
           arch.getNAodIntermediateLevels();
  }
  if (aodMovesNeighbor.empty()) {
    return getMaxOffset(dim, x, sign) < arch.getNAodIntermediateLevels();
  }
  return getMaxOffset(dim, x, sign) + getMaxOffset(dim, neighborX, sign) <
         arch.getNAodIntermediateLevels();
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
    if (std::find(qubitsMove.begin(), qubitsMove.end(), move.first) ==
        qubitsMove.end()) {
      qubitsMove.push_back(move.first);
    }
    if (std::find(qubitsMove.begin(), qubitsMove.end(), move.second) ==
        qubitsMove.end()) {
      qubitsMove.push_back(move.second);
    }
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
