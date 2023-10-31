//
// Created by Ludwig Schmid on 05.10.23.
//

#include "namap/NeutralAtomMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "namap/NeutralAtomUtils.hpp"
#include "utils.hpp"

namespace qc {
QuantumComputation qc::NeutralAtomMapper::map(qc::QuantumComputation& qc) {
  qc::CircuitOptimizer::removeFinalMeasurements(qc);
  this->dag     = qc::CircuitOptimizer::constructDAG(qc);
  this->mapping = Mapping(qc.getNqubits(), InitialMapping::Identity);
  for (auto& i : dag) {
    this->frontLayerIterators.push_back(i.begin());
    this->lookaheadOffsets.push_back(0);
    this->lookaheadCandidates.emplace_back();
    this->frontCandidates.emplace_back();
  }

  if (dag.size() > arch.getNqubits()) {
    throw std::runtime_error("More qubits in circuit than in architecture");
  }

  std::cout << "test" << std::endl;

  this->parameters.lookaheadWeightSwaps = 0.0;
  this->parameters.lookaheadWeightMoves = 0.1;
  this->parameters.decay                = 1;
  this->parameters.shuttlingTimeWeight  = 0;
  //  this->parameters.shuttlingMakeExecutableBonus = arch.getNcolumns();
  this->parameters.shuttlingMakeExecutableBonus = 1;
  this->verbose                                 = false;

  //   precompute all distances

  this->decayWeights.reserve(this->arch.getNcolumns());
  for (uint32_t i = this->arch.getNcolumns(); i > 0; --i) {
    this->decayWeights.push_back(std::exp(-this->parameters.decay * i));
  }

  createFrontLayer();
  //  printLayers();
  updateLookaheadLayerByQubit();
  if (this->verbose) {
    printLayers();
  }

  //  hardwareQubits.move(0, 15, arch);
  // probably need
  // updateFrontLayerByQubit(list with all qubits) to trigger commutation

  // TODO: find SWAP
  // TODO: get list of front gates that can be executed now
  auto i = 0;
  while (!this->frontLayer.empty()) {
    i++;
    if (this->verbose) {
      std::cout << "iteration " << i << std::endl;
    }
    std::set<std::unique_ptr<Operation>*> gatesToExecute;
    while (gatesToExecute.empty()) {
      //      auto bestSwap = findBestSwap();
      //      updateMapping(bestSwap);
      auto bestMove = findBestAtomMove();
      updateMappingMove(bestMove);
      gatesToExecute = getExecutableGates();
    }
    updateFrontLayerByGate(gatesToExecute);
    //    printLayers();
    updateLookaheadLayerByQubit();
    if (this->verbose) {
      printLayers();
    }
  }
  std::cout << "nSwaps: " << nSwaps << std::endl;
  std::cout << "nMoves: " << nMoves << std::endl;
  qc::CircuitOptimizer::decomposeSWAP(this->mappedQc, false);
  AodScheduler scheduler(this->arch);
  auto         scheduledQc = scheduler.schedule(this->mappedQc);
  return scheduledQc;
}

void qc::NeutralAtomMapper::createFrontLayer() {
  // vector indicating if qubit is free
  std::vector<Qubit> allQubits;
  for (uint32_t i = 0; i < this->dag.size(); i++) {
    allQubits.push_back(i);
    this->frontQubitsToUpdate.insert(i);
    this->lookaheadQubitsToUpdate.insert(i);
  }
  updateFrontLayerByQubit();
  updateLookaheadLayerByQubit();
  lookaheadQubitsToUpdate.clear();
}

void qc::NeutralAtomMapper::updateFrontLayerByGate(
    std::set<std::unique_ptr<Operation>*>& gatesToExecute) {
  for (auto* gate : gatesToExecute) {
    // add other gate qubits to be checked
    for (auto qubit : (*gate)->getUsedQubits()) {
      this->frontQubitsToUpdate.insert(qubit);
      this->lookaheadOffsets[qubit]--;
      this->frontLayerIterators[qubit]++;
    }
    // remove from current FrontLayer
    mapGate(gate);
    // update DAG iterators
    // remove from FrontLayer
    this->frontLayer.erase(
        std::find(this->frontLayer.begin(), this->frontLayer.end(), gate));
  }
  updateFrontLayerByQubit();
}

void qc::NeutralAtomMapper::updateFrontLayerByQubit() {
  findFrontCandidates();
  updateFrontLayerByCandidates();
}

void qc::NeutralAtomMapper::findFrontCandidates() {
  for (auto& qubit : this->frontQubitsToUpdate) {
    auto tempIter = this->frontLayerIterators[qubit];
    while (tempIter != this->dag[qubit].end()) {
      auto* opPointer = *tempIter;
      if (opPointer->get()->getUsedQubits().size() == 1) {
        mapGate(opPointer);
        this->frontLayerIterators[qubit]++;
        tempIter++;
      } else {
        // continue if following gates commute
        bool commutes = true;
        while (commutes && tempIter != this->dag[qubit].end()) {
          auto* nextOpPointer = *tempIter;
          commutes            = commutesWith(this->frontLayer, nextOpPointer) &&
                     commutesWith(this->frontCandidates[qubit], nextOpPointer);
          if (commutes) {
            if (nextOpPointer->get()->getUsedQubits().size() == 1) {
              mapGate(nextOpPointer);
            } else { // not executable but commutes
              this->frontCandidates[qubit].insert(nextOpPointer);
            }
            tempIter++;
          }
        }
        break;
      }
    }
  }
}

void qc::NeutralAtomMapper::updateFrontLayerByCandidates() {
  // create deep copy of frontQubitsToUpdate
  const std::vector<Qubit> tempQubitsToUpdate(this->frontQubitsToUpdate.begin(),
                                              this->frontQubitsToUpdate.end());
  this->frontQubitsToUpdate.clear();
  for (auto& qubit : tempQubitsToUpdate) {
    std::vector<std::unique_ptr<qc::Operation>*> toRemove;
    for (auto* opPointer : this->frontCandidates[qubit]) {
      bool  inFrontLayer = true;
      auto* op           = opPointer->get();
      if (op->getType() == qc::OpType::I) {
        continue;
      }
      for (const auto& opQubit : op->getUsedQubits()) {
        if (qubit == opQubit) {
          continue;
        }
        if (std::find(this->frontCandidates[opQubit].begin(),
                      this->frontCandidates[opQubit].end(),
                      opPointer) == this->frontCandidates[opQubit].end()) {
          inFrontLayer = false;
          break;
        }
      }
      if (inFrontLayer) {
        // if executable execute directly
        if (isExecutable(opPointer)) {
          mapGate(opPointer);
          for (const auto& opQubit : op->getUsedQubits()) {
            this->frontLayerIterators[opQubit]++;
            this->frontQubitsToUpdate.insert(opQubit);
          }
          // remove from lookahead layer if there
          if (this->lookaheadLayer.find(opPointer) !=
              this->lookaheadLayer.end()) {
            this->lookaheadLayer.erase(opPointer);
          }
        } else {
          addToFrontLayer(opPointer);
          // move lookahead
          for (const auto& opQubit : op->getUsedQubits()) {
            this->lookaheadOffsets[opQubit]++;
          }
        }
        // remove from candidacy of other qubits
        for (const auto& opQubit : op->getUsedQubits()) {
          this->lookaheadQubitsToUpdate.insert(opQubit);
          if (qubit == opQubit) {
            continue;
          }
          this->frontCandidates[opQubit].erase(opPointer);
        }

        // save to remove from candidacy of this qubit
        toRemove.push_back(opPointer);
      }
    }
    // remove from candidacy of this qubit
    // has to be done now to not change iterating list
    for (auto* opPointer : toRemove) {
      this->frontCandidates[qubit].erase(opPointer);
    }
  }
  // if gates have been executed directly refresh front layer
  if (!this->frontQubitsToUpdate.empty()) {
    updateFrontLayerByQubit();
  }
}

bool qc::NeutralAtomMapper::commutesWith(
    GateList& layer, std::unique_ptr<qc::Operation>* opPointer) {
  return std::all_of(layer.begin(), layer.end(),
                     [&opPointer](const auto& frontOpPointer) {
                       return commute(opPointer, frontOpPointer);
                     });
}

bool qc::NeutralAtomMapper::commute(
    std::unique_ptr<qc::Operation>* opPointer1,
    std::unique_ptr<qc::Operation>* opPointer2) {
  auto* op1 = opPointer1->get();
  auto* op2 = opPointer2->get();

  // gates acting of disjoint sets of qubits commute
  auto usedQubits1 = op1->getUsedQubits();
  auto usedQubits2 = op2->getUsedQubits();
  if (std::all_of(usedQubits1.begin(), usedQubits1.end(),
                  [&usedQubits2](const auto& qubit) {
                    return std::find(usedQubits2.begin(), usedQubits2.end(),
                                     qubit) == usedQubits2.end();
                  })) {
    return true;
  }

  // single qubit gates commute
  if (op1->getUsedQubits().size() == 1 && op2->getUsedQubits().size() == 1) {
    return true;
  }

  // if targets are same, it commutes for same type
  if (op1->getTargets() == op2->getTargets() &&
      op1->getType() == op2->getType()) {
    return true;
  }

  // if controls act on same qubits, it commutes
  auto controls1 = op1->getControls();
  auto controls2 = op2->getControls();
  if (!controls1.empty() && !controls2.empty()) {
    if (std::all_of(controls1.begin(), controls1.end(),
                    [&controls2](const auto& control1) {
                      return std::any_of(controls2.begin(), controls2.end(),
                                         [&control1](const auto& control2) {
                                           return control1.qubit ==
                                                  control2.qubit;
                                         });
                    })) {
      return true;
    }
  }

  // control and single qubit z-rz commutes

  if (commuteSingleQubitZAndControl(opPointer1, opPointer2) ||
      commuteSingleQubitZAndControl(opPointer2, opPointer1)) {
    return true;
  }
  return false;
}

bool qc::NeutralAtomMapper::commuteSingleQubitZAndControl(
    std::unique_ptr<qc::Operation>* opPointer1,
    std::unique_ptr<qc::Operation>* opPointer2) {
  auto* op1 = opPointer1->get();
  auto* op2 = opPointer2->get();
  if (op1->getType() == qc::OpType::Z || op1->getType() == qc::OpType::RZ) {
    if (op1->getControls().empty()) {
      auto op1Target   = op1->getTargets().front();
      auto op2Controls = op2->getControls();
      if (std::any_of(op2Controls.begin(), op2Controls.end(),
                      [&op1Target](const auto& control) {
                        return op1Target == control.qubit;
                      })) {
        return true;
      }
    }
  }
  return false;
}

void qc::NeutralAtomMapper::updateLookaheadLayerByQubit() {
  // find all multi qubit gates within the lookahead depth

  findLookaheadCandidates();
  updateLookaheadLayerByCandidates();

  this->lookaheadQubitsToUpdate.clear();
}

void qc::NeutralAtomMapper::updateLookaheadLayerByCandidates() {
  for (auto& qubit : this->lookaheadQubitsToUpdate) {
    std::vector<std::unique_ptr<qc::Operation>*> toRemove;
    for (auto* opPointer : this->lookaheadCandidates[qubit]) {
      bool  inLookahead = true;
      auto* op          = opPointer->get();
      for (const auto& opQubit : op->getUsedQubits()) {
        if (qubit == opQubit) {
          continue;
        }
        if (std::find(this->lookaheadCandidates[opQubit].begin(),
                      this->lookaheadCandidates[opQubit].end(),
                      opPointer) == this->lookaheadCandidates[opQubit].end()) {
          inLookahead = false;
          break;
        }
      }
      if (inLookahead) {
        this->lookaheadLayer.insert(opPointer);
        // remove from candidacy of other qubits
        for (const auto& opQubit : op->getUsedQubits()) {
          if (qubit == opQubit) {
            continue;
          }
          this->lookaheadCandidates[opQubit].erase(opPointer);
        }
        // save to remove from candidacy of this qubit
        toRemove.push_back(opPointer);
      }
    }
    // remove from candidacy of this qubit
    for (auto* opPointer : toRemove) {
      this->lookaheadCandidates[qubit].erase(opPointer);
    }
  }
}

void qc::NeutralAtomMapper::findLookaheadCandidates() {
  // find all multi qubit gates within the lookahead depth
  // for each qubit that needs to be updated
  for (auto& qubit : this->lookaheadQubitsToUpdate) {
    // init Iterator
    if (this->frontLayerIterators[qubit] + this->lookaheadOffsets[qubit] >=
        this->dag[qubit].end()) {
      continue;
    }
    auto lookaheadIter =
        this->frontLayerIterators[qubit] + this->lookaheadOffsets[qubit];

    // iterate until end or lookahead depth reached
    uint32_t foundMultiQubitGate = 0;
    while (lookaheadIter != this->dag[qubit].end() &&
           foundMultiQubitGate < this->lookaheadDepth) {
      auto* opPointer = *lookaheadIter;
      auto* op        = opPointer->get();
      if (op->getType() == qc::OpType::I) {
        lookaheadIter++;
        continue;
      }
      if (op->getUsedQubits().size() > 1) {
        foundMultiQubitGate++;
        this->lookaheadCandidates[qubit].insert(opPointer);
        lookaheadIter++;
        // continue if following gates commute
        // only check last gate (only approximate)
        bool commutes = true;
        while (commutes && lookaheadIter != this->dag[qubit].end()) {
          auto* nextOpPointer = *lookaheadIter;
          auto* nextOp        = nextOpPointer->get();
          commutes            = commute(opPointer, nextOpPointer);
          if (commutes) {
            lookaheadIter++;
            if (nextOp->getUsedQubits().size() > 1) {
              this->lookaheadCandidates[qubit].insert(nextOpPointer);
            }
          } else {
            commutes = false;
          }
        }
      }
      lookaheadIter++;
    }
  }
}

void qc::NeutralAtomMapper::mapGate(std::unique_ptr<qc::Operation>* op) {
  if (op->get()->getType() == qc::OpType::I) {
    return;
  }
  // convert circuit qubits to CoordIndex and append to mappedQc
  this->mapping.mapToHwQubits(op);
  this->hardwareQubits.mapToCoordIdx(op);
  this->mappedQc.emplace_back((*op)->clone());
  if (this->verbose) {
    std::cout << "mapped " << (*op)->getName() << " ";
    for (auto qubit : (*op)->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << "\n";
  }

  op->get()->setGate(qc::OpType::I);
}

bool qc::NeutralAtomMapper::isExecutable(
    std::unique_ptr<qc::Operation>* opPointer) {
  auto* op          = opPointer->get();
  auto  usedQubits  = op->getUsedQubits();
  auto  nUsedQubits = usedQubits.size();
  if (nUsedQubits == 1) {
    return true;
  }
  std::set<Qubit> usedHwQubits;
  for (auto qubit : usedQubits) {
    usedHwQubits.insert(this->mapping.getHwQubit(qubit));
  }
  return this->hardwareQubits.getTotalDistance(usedHwQubits) == 0;
}

void qc::NeutralAtomMapper::addToFrontLayer(
    std::unique_ptr<qc::Operation>* opPointer) {
  this->frontLayer.insert(opPointer);
  // remove from lookahead layer if there
  if (this->lookaheadLayer.find(opPointer) != this->lookaheadLayer.end()) {
    this->lookaheadLayer.erase(opPointer);
  }
}

void qc::NeutralAtomMapper::printLayers() {
  std::cout << "f: ";
  for (auto* opPointer : this->frontLayer) {
    auto* op = opPointer->get();
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "l: ";
  for (auto* opPointer : this->lookaheadLayer) {
    auto* op = opPointer->get();
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

GateList NeutralAtomMapper::getExecutableGates() {
  GateList executableGates;
  for (auto* opPointer : this->frontLayer) {
    if (isExecutable(opPointer)) {
      executableGates.insert(opPointer);
    }
  }
  return executableGates;
}

void NeutralAtomMapper::updateMapping(qc::Swap swap) {
  nSwaps++;
  // save to lastSwaps
  this->lastBlockedQubits.push_back(
      this->hardwareQubits.getBlockedQubits({swap.first, swap.second}));
  if (this->lastBlockedQubits.size() > this->arch.getNcolumns()) {
    this->lastBlockedQubits.pop_front();
  }
  this->mapping.swap(swap);
  // convert circuit qubits to CoordIndex and append to mappedQc
  auto idxFirst =
      this->hardwareQubits.getCoordIndex(this->mapping.getHwQubit(swap.first));
  auto idxSecond =
      this->hardwareQubits.getCoordIndex(this->mapping.getHwQubit(swap.second));
  this->mappedQc.swap(idxFirst, idxSecond);
  if (this->verbose) {
    std::cout << "swapped " << swap.first << " " << swap.second;
    std::cout << "  logical qubits: ";
    if (this->mapping.isMapped(swap.first)) {
      std::cout << this->mapping.getCircQubit(swap.first);
    } else {
      std::cout << "not mapped";
    }
    if (this->mapping.isMapped(swap.second)) {
      std::cout << " " << this->mapping.getCircQubit(swap.second);
    } else {
      std::cout << " not mapped";
    }
    std::cout << std::endl;
  }
}

void NeutralAtomMapper::updateMappingMove(qc::AtomMove move) {
  this->lastMoves.push_back(move);
  if (this->lastMoves.size() > this->arch.getNcolumns()) {
    this->lastMoves.pop_front();
  }
  mappedQc.move(move.first, move.second);
  this->mappedQc.move(move.first, move.second);
  auto toMoveHwQubit = this->hardwareQubits.getHwQubit(move.first);
  this->hardwareQubits.move(toMoveHwQubit, move.second);
  if (verbose) {
    std::cout << "moved " << move.first << " to " << move.second << " HwQubit "
              << toMoveHwQubit << " CircQubit "
              << this->mapping.getCircQubit(toMoveHwQubit) << std::endl;
  }
  nMoves++;
}

qc::Swap qc::NeutralAtomMapper::findBestSwap() {
  auto swaps = getAllPossibleSwaps();
  // evaluate swaps based on cost function
  std::vector<std::pair<Swap, fp>> swapCosts;
  swapCosts.reserve(swaps.size());
  for (const auto& swap : swaps) {
    swapCosts.emplace_back(swap, distanceCost(swap));
  }
  // get swap of minimal cost
  auto bestSwap = std::min_element(swapCosts.begin(), swapCosts.end(),
                                   [](const auto& swap1, const auto& swap2) {
                                     return swap1.second < swap2.second;
                                   });
  // none of the swaps helps
  if (bestSwap->second == 0) {
    throw std::runtime_error("Qubit not reachable -> no swap helps");
  }
  return bestSwap->first;
}

std::set<qc::Swap> qc::NeutralAtomMapper::getAllPossibleSwaps() {
  std::set<Swap> swaps;
  for (const auto& op : this->frontLayer) {
    for (const auto& qubit : op->get()->getUsedQubits()) {
      auto nearbySwaps =
          this->hardwareQubits.getNearbySwaps(this->mapping.getHwQubit(qubit));
      for (const auto& swap : nearbySwaps) {
        swaps.insert(swap);
      }
    }
  }
  return swaps;
}

fp NeutralAtomMapper::distanceCost(const qc::Swap& swap) {
  // compute the change in total distance
  auto distanceChangeFront = distancePerLayer(swap, this->frontLayer) /
                             static_cast<fp>(this->frontLayer.size());
  auto distanceChangeLookahead = distancePerLayer(swap, this->lookaheadLayer) /
                                 static_cast<fp>(this->lookaheadLayer.size());
  auto cost = parameters.lookaheadWeightSwaps * distanceChangeLookahead +
              distanceChangeFront;
  //  compute the last time one of the swap qubits was used
  if (this->parameters.decay != 0) {
    uint32_t idxLastUsed = 0;
    for (uint32_t i = 0; i < this->lastBlockedQubits.size(); ++i) {
      if (this->lastBlockedQubits[i].find(swap.first) !=
              this->lastBlockedQubits[i].end() ||
          this->lastBlockedQubits[i].find(swap.second) !=
              this->lastBlockedQubits[i].end()) {
        idxLastUsed = i;
        break;
      }
    }
    cost *= this->decayWeights[idxLastUsed];
  }
  return cost;
}

fp NeutralAtomMapper::distancePerLayer(const qc::Swap& swap, GateList& layer) {
  fp distBefore = 0;
  fp distAfter  = 0;
  fp distChange = 0;
  // bring close only ontil swap distance =0, bring exact to the exact position
  std::vector<std::pair<HwQubit, HwQubit>> swapCloseBy;
  std::vector<std::pair<HwQubit, HwQubit>> moveExact;
  for (const auto& gate : layer) {
    auto usedQubits = gate->get()->getUsedQubits();
    if (usedQubits.size() == 2) {
      swapCloseBy.emplace_back(this->mapping.getHwQubit(*usedQubits.begin()),
                               this->mapping.getHwQubit(*usedQubits.rbegin()));
    } else {
      // TODO get list of origin -> target from optimal positions for multi
      // qubit gates
      throw std::runtime_error(
          "Multi qubit gates distance not yet implemented");
    }
  }
  // bring qubits together to execute gate
  for (const auto& [q1, q2] : swapCloseBy) {
    // distance before
    distBefore = this->hardwareQubits.getSwapDistance(q1, q2);
    if (distBefore == std::numeric_limits<fp>::infinity()) {
      continue;
    }
    // do swap
    if (q1 == swap.first) {
      distAfter = this->hardwareQubits.getSwapDistance(swap.second, q2);
    } else if (q2 == swap.second) {
      distAfter = this->hardwareQubits.getSwapDistance(q1, swap.first);
    } else if (q1 == swap.second) {
      distAfter = this->hardwareQubits.getSwapDistance(swap.first, q2);
    } else if (q2 == swap.first) {
      distAfter = this->hardwareQubits.getSwapDistance(q1, swap.second);
    } else {
      continue;
    }
    distChange += distAfter - distBefore;
  }

  // move qubits to the exact position for multi-qubit gates
  for (const auto& [origin, destination] : moveExact) {
    distBefore =
        this->hardwareQubits.getSwapDistance(origin, destination, false);
    if (distBefore == std::numeric_limits<fp>::infinity()) {
      continue;
    }
    if (origin == swap.first) {
      distAfter =
          this->hardwareQubits.getSwapDistance(swap.second, destination, false);
    } else if (origin == swap.second) {
      distAfter =
          this->hardwareQubits.getSwapDistance(swap.first, destination, false);
    }
    distChange += distAfter - distBefore;
  }

  return distChange;
}

AtomMove NeutralAtomMapper::findBestAtomMove() {
  auto moveCombs = getAllPossibleMoveCombinations();

  std::vector<std::pair<AtomMove, fp>> moveCosts;
  moveCosts.reserve(moveCombs.size());
  for (const auto& moveComb : moveCombs) {
    moveCosts.emplace_back(*moveComb.rbegin(), moveCostComb(moveComb));
  }

  // get move of minimal cost
  auto bestMove = std::min_element(moveCosts.begin(), moveCosts.end(),
                                   [](const auto& move1, const auto& move2) {
                                     return move1.second < move2.second;
                                   });
  return bestMove->first;
}

fp NeutralAtomMapper::moveCostComb(const qc::MoveComb& moveComb) {
  fp costComb = 0;
  for (const auto& move : moveComb) {
    costComb += moveCost(move);
  }
  return costComb;
}

fp NeutralAtomMapper::moveCost(const AtomMove& move) {
  fp   cost      = 0;
  auto frontCost = moveDistancePerLayer(move, this->frontLayer) /
                   static_cast<fp>(this->frontLayer.size());
  cost += frontCost;
  if (!lookaheadLayer.empty()) {
    auto lookaheadCost = moveDistancePerLayer(move, this->lookaheadLayer) /
                         static_cast<fp>(this->lookaheadLayer.size());
    cost += parameters.lookaheadWeightMoves * lookaheadCost;
  }
  auto parallelCost = parameters.shuttlingTimeWeight * parallelMoveCost(move);
  cost += parallelCost;
  return cost;
}

fp NeutralAtomMapper::moveDistancePerLayer(const qc::AtomMove& move,
                                           qc::GateList&       layer) {
  // compute cost assuming the move was applied
  fp   distChange    = 0;
  auto toMoveHwQubit = this->hardwareQubits.getHwQubit(move.first);
  if (this->mapping.isMapped(toMoveHwQubit)) {
    auto toMoveCircuitQubit = this->mapping.getCircQubit(toMoveHwQubit);
    for (const auto& gate : layer) {
      auto usedQubits = gate->get()->getUsedQubits();
      if (usedQubits.find(toMoveCircuitQubit) != usedQubits.end()) {
        // check distance reduction
        fp distanceBefore = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          distanceBefore += this->hardwareQubits.getSwapDistance(toMoveHwQubit,
                                                                 hwQubit, true);
        }
        fp distanceAfter = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          distanceAfter +=
              this->hardwareQubits.getSwapDistanceMove(move.second, hwQubit);
        }
        distChange += distanceAfter - distanceBefore;

        // add bonus if distance is 0 (not requiring another shuttling)
        if (distanceAfter == 0) {
          distChange -= parameters.shuttlingMakeExecutableBonus;
        }
      }
    }
  }
  return distChange;
}

fp NeutralAtomMapper::parallelMoveCost(const qc::AtomMove& move) {
  fp   parallelCost = 0;
  auto moveVector   = this->arch.getVector(move.first, move.second);
  std::vector<CoordIndex> lastEndingCoords;
  if (this->lastMoves.empty()) {
    parallelCost += arch.getVectorShuttlingTime(moveVector);
  }
  for (const auto& lastMove : this->lastMoves) {
    lastEndingCoords.push_back(lastMove.second);
    // decide of shuttling can be done in parallel
    auto lastMoveVector = this->arch.getVector(lastMove.first, lastMove.second);
    if (moveVector.overlap(lastMoveVector)) {
      if (moveVector.direction != lastMoveVector.direction) {
        parallelCost += arch.getVectorShuttlingTime(moveVector);
      } else {
        // check if move can be done in parallel
        if (moveVector.include(lastMoveVector)) {
          parallelCost += arch.getVectorShuttlingTime(moveVector);
        }
      }
    }
  }
  // check if move can use AOD atom from last moves
  if (std::find(lastEndingCoords.begin(), lastEndingCoords.end(), move.first) ==
      lastEndingCoords.end()) {
    parallelCost += arch.getShuttlingTime(OpType::AodActivate) +
                    arch.getShuttlingTime(OpType::AodDeactivate);
  }
  return parallelCost;
}

std::set<MoveComb> NeutralAtomMapper::getAllPossibleMoveCombinations() {
  std::set<MoveComb> moves;
  for (const auto& op : this->frontLayer) {
    auto usedQubits   = op->get()->getUsedQubits();
    auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
    if (usedHwQubits.size() == 2) {
      auto q1      = *usedHwQubits.begin();
      auto q2      = *usedHwQubits.rbegin();
      auto moves12 = getNearbyMoveCombinations(q1, q2);
      auto moves21 = getNearbyMoveCombinations(q2, q1);
      moves.insert(moves12.begin(), moves12.end());
      moves.insert(moves21.begin(), moves21.end());
      // TODO test if one should always include move away in both directions
      if (moves.empty()) {
        auto moves1 = getMoveAwayCombinations(q1, q2);
        auto moves2 = getMoveAwayCombinations(q2, q1);
        moves.insert(moves1.begin(), moves1.end());
        moves.insert(moves2.begin(), moves2.end());
      }
    } else {
      throw std::runtime_error(
          "Multi qubit gates distance not yet implemented");
    }
  }
  return moves;
}

std::set<MoveComb>
NeutralAtomMapper::getNearbyMoveCombinations(HwQubit start, HwQubit target) {
  // Finds free coords near target and moves start to them
  auto               startCoord = this->hardwareQubits.getCoordIndex(start);
  std::set<MoveComb> moveCombinations;
  auto nearbyFreeCoords = this->hardwareQubits.getNearbyFreeCoordinates(target);
  for (const auto& coord : nearbyFreeCoords) {
    const AtomMove move = {startCoord, coord};
    moveCombinations.insert({move});
  }
  return moveCombinations;
}

std::set<MoveComb>
NeutralAtomMapper::getMoveAwayCombinations(qc::HwQubit start,
                                           qc::HwQubit target) {
  // Finds all possibilities to move a nearby atom away from target
  // and then move start to the now free coord
  std::set<MoveComb> moveCombinations;

  auto startCoord         = this->hardwareQubits.getCoordIndex(start);
  auto nearbyTargetQubits = this->hardwareQubits.getNearbyQubits(target);
  for (const auto& moveAwayQubit : nearbyTargetQubits) {
    auto moveAwayStart = this->hardwareQubits.getCoordIndex(moveAwayQubit);
    auto const originalVector = this->arch.getVector(startCoord, moveAwayStart);
    auto const originalDirection = originalVector.direction;
    auto       moveAwayTargets   = this->hardwareQubits.findClosestFreeCoord(
        moveAwayQubit, originalDirection);
    for (const auto& moveAwayTarget : moveAwayTargets) {
      const AtomMove move     = {startCoord, moveAwayStart};
      const AtomMove moveAway = {moveAwayStart, moveAwayTarget};
      moveCombinations.insert({move, moveAway});
    }
  }
  return moveCombinations;
}

} // namespace qc
