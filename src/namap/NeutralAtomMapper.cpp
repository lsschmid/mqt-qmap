//
// Created by Ludwig Schmid on 05.10.23.
//

#include "namap/NeutralAtomMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "namap/NeutralAtomUtils.hpp"
#include "utils.hpp"

namespace qc {
QuantumComputation qc::NeutralAtomMapper::map(qc::QuantumComputation& qc,
                                              InitialMapping initialMapping,
                                              bool           verboseArg) {
  qc::CircuitOptimizer::removeFinalMeasurements(qc);
  this->dag     = qc::CircuitOptimizer::constructDAG(qc);
  this->mapping = Mapping(qc.getNqubits(), initialMapping);
  for (auto& i : dag) {
    this->frontLayerIterators.push_back(i.begin());
    this->lookaheadOffsets.push_back(0);
    this->lookaheadCandidates.emplace_back();
    this->frontCandidates.emplace_back();
  }

  if (dag.size() > arch.getNqubits()) {
    throw std::runtime_error("More qubits in circuit than in architecture");
  }

  std::cout << "test" << '\n';

  this->verbose = verboseArg;

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

  // TODO: find SWAP
  // TODO: get list of front gates that can be executed now
  auto i = 0;
  while (!this->frontLayerGate.empty() || !this->frontLayerShuttling.empty()) {
    // first do all gate based Mapping
    while (!this->frontLayerGate.empty()) {
      std::set<std::unique_ptr<Operation>*> gatesToExecute;
      Swap lastSwap = {std::numeric_limits<uint32_t>::max(),
                       std::numeric_limits<uint32_t>::max()};
      while (gatesToExecute.empty()) {
        ++i;
        if (this->verbose) {
          std::cout << "iteration " << i << '\n';
        }
        auto bestSwap = findBestSwap();
        if (bestSwap.first == std::numeric_limits<uint32_t>::max()) {
          break;
        }
        lastSwap = bestSwap;
        updateMapping(bestSwap);
        //      auto bestMove = findBestAtomMove();
        //      updateMappingMove(bestMove);
        gatesToExecute = getExecutableGates();
      }
      updateFrontLayerByGate(gatesToExecute);
      updateLookaheadLayerByQubit();
      if (this->verbose) {
        printLayers();
      }
    }
    while (!this->frontLayerShuttling.empty()) {
      std::set<std::unique_ptr<Operation>*> gatesToExecute;
      while (gatesToExecute.empty()) {
        ++i;
        if (this->verbose) {
          std::cout << "iteration " << i << '\n';
        }
        auto bestMove = findBestAtomMove();
        updateMappingMove(bestMove);
        gatesToExecute = getExecutableGates();
      }
      updateFrontLayerByGate(gatesToExecute);
      updateLookaheadLayerByQubit();
      if (this->verbose) {
        printLayers();
      }
    }
  }
  std::cout << "nSwaps: " << nSwaps << '\n';
  std::cout << "nMoves: " << nMoves << '\n';
  return this->mappedQc;
}

QuantumComputation NeutralAtomMapper::mapAod(qc::QuantumComputation& qc) {
  qc::CircuitOptimizer::decomposeSWAP(this->mappedQc, false);
  AodScheduler scheduler(this->arch);
  return scheduler.schedule(qc);
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
    if (this->frontLayerGate.find(gate) != this->frontLayerGate.end()) {
      this->frontLayerGate.erase(std::find(this->frontLayerGate.begin(),
                                           this->frontLayerGate.end(), gate));
    }
    if (this->frontLayerShuttling.find(gate) !=
        this->frontLayerShuttling.end()) {
      this->frontLayerShuttling.erase(
          std::find(this->frontLayerShuttling.begin(),
                    this->frontLayerShuttling.end(), gate));
    }
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
          commutes =
              commutesWithAtQubit(this->frontLayerGate, nextOpPointer, qubit) &&
              commutesWithAtQubit(this->frontLayerShuttling, nextOpPointer,
                                  qubit) &&
              commutesWithAtQubit(this->frontCandidates[qubit], nextOpPointer,
                                  qubit);
          if (commutes) {
            if (nextOpPointer->get()->getUsedQubits().size() == 1) {
              mapGate(nextOpPointer);
            } else { // not executable but commutes
              if (this->frontLayerGate.find(nextOpPointer) ==
                      this->frontLayerGate.end() &&
                  this->frontLayerShuttling.find(nextOpPointer) ==
                      this->frontLayerShuttling.end()) {
                this->frontCandidates[qubit].insert(nextOpPointer);
              }
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
        toRemove.push_back(opPointer);
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
          if (this->lookaheadLayerGate.find(opPointer) !=
              this->lookaheadLayerGate.end()) {
            this->lookaheadLayerGate.erase(opPointer);
          }
          if (this->lookaheadLayerShuttling.find(opPointer) !=
              this->lookaheadLayerShuttling.end()) {
            this->lookaheadLayerShuttling.erase(opPointer);
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

bool qc::NeutralAtomMapper::commutesWithAtQubit(
    const qc::GateList& layer, const std::unique_ptr<qc::Operation>* opPointer,
    const Qubit& qubit) {
  return std::all_of(layer.begin(), layer.end(),
                     [&opPointer, &qubit](const auto& frontOpPointer) {
                       return commuteAtQubit(opPointer, frontOpPointer, qubit);
                     });
}

bool qc::NeutralAtomMapper::commuteAtQubit(
    const std::unique_ptr<qc::Operation>* opPointer1,
    const std::unique_ptr<qc::Operation>* opPointer2, const qc::Qubit& qubit) {
  auto* op1 = opPointer1->get();
  auto* op2 = opPointer2->get();

  // single qubit gates commute
  if (op1->getUsedQubits().size() == 1 && op2->getUsedQubits().size() == 1) {
    return true;
  }

  if (op1->getType() == qc::OpType::I || op2->getType() == qc::OpType::I) {
    return true;
  }

  // commutes at qubit if at least one of the two gates does not use qubit
  auto usedQubits1 = op1->getUsedQubits();
  auto usedQubits2 = op2->getUsedQubits();
  if (usedQubits1.find(qubit) == usedQubits1.end() ||
      usedQubits2.find(qubit) == usedQubits2.end()) {
    return true;
  }

  // for two-qubit gates, check if they commute at qubit
  // commute if both are controlled at qubit or operation on qubit is same
  // check controlles
  if (op1->getControls().find(qubit) != op1->getControls().end() &&
      op2->getControls().find(qubit) != op2->getControls().end()) {
    return true;
  }
  // check targets
  if (std::find(op1->getTargets().begin(), op1->getTargets().end(), qubit) !=
          op1->getTargets().end() &&
      (std::find(op2->getTargets().begin(), op2->getTargets().end(), qubit) !=
       op2->getTargets().end()) &&
      (op1->getType() == op2->getType())) {
    return true;
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
  for (const auto& qubit : this->lookaheadQubitsToUpdate) {
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
        addToLookaheadLayer(opPointer);
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
          commutes            = commuteAtQubit(opPointer, nextOpPointer, qubit);
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
  if (this->verbose) {
    std::cout << "mapped " << (*op)->getName() << " ";
    for (auto qubit : (*op)->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << "\n";
  }
  // convert circuit qubits to CoordIndex and append to mappedQc
  auto opClone = (*op)->clone();
  this->mapping.mapToHwQubits(&opClone);
  this->hardwareQubits.mapToCoordIdx(&opClone);
  this->mappedQc.emplace_back(opClone);

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
  if (swapGateBetter(opPointer)) {
    this->frontLayerGate.insert(opPointer);
    // remove from lookahead layer if there
    if (this->lookaheadLayerGate.find(opPointer) !=
        this->lookaheadLayerGate.end()) {
      this->lookaheadLayerGate.erase(opPointer);
    }
  } else {
    this->frontLayerShuttling.insert(opPointer);
    // remove from lookahead layer if there
    if (this->lookaheadLayerShuttling.find(opPointer) !=
        this->lookaheadLayerShuttling.end()) {
      this->lookaheadLayerShuttling.erase(opPointer);
    }
  }
  // remove from front candidates not allowed here
  // to prevent modify while iterating
}

void NeutralAtomMapper::addToLookaheadLayer(
    std::unique_ptr<qc::Operation>* opPointer) {
  if (swapGateBetter(opPointer)) {
    this->lookaheadLayerGate.insert(opPointer);
  } else {
    this->lookaheadLayerShuttling.insert(opPointer);
  }
}

void qc::NeutralAtomMapper::printLayers() {
  std::cout << "f,g: ";
  for (auto* opPointer : this->frontLayerGate) {
    auto* op = opPointer->get();
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << "f,s: ";
  for (auto* opPointer : this->frontLayerShuttling) {
    auto* op = opPointer->get();
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << "l,g: ";
  for (auto* opPointer : this->lookaheadLayerGate) {
    auto* op = opPointer->get();
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
  std::cout << "l,g: ";
  for (auto* opPointer : this->lookaheadLayerShuttling) {
    auto* op = opPointer->get();
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

GateList NeutralAtomMapper::getExecutableGates() {
  GateList executableGates;
  for (auto* opPointer : this->frontLayerGate) {
    if (isExecutable(opPointer)) {
      executableGates.insert(opPointer);
    }
  }
  for (auto* opPointer : this->frontLayerShuttling) {
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
  auto idxFirst  = this->hardwareQubits.getCoordIndex(swap.first);
  auto idxSecond = this->hardwareQubits.getCoordIndex(swap.second);
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
    std::cout << '\n';
  }
}

void NeutralAtomMapper::updateMappingMove(qc::AtomMove move) {
  this->lastMoves.push_back(move);
  if (this->lastMoves.size() > this->arch.getNcolumns()) {
    this->lastMoves.pop_front();
  }
  mappedQc.move(move.first, move.second);
  auto toMoveHwQubit = this->hardwareQubits.getHwQubit(move.first);
  this->hardwareQubits.move(toMoveHwQubit, move.second);
  if (verbose) {
    std::cout << "moved " << move.first << " to " << move.second;
    if (this->mapping.isMapped(toMoveHwQubit)) {
      std::cout << "  logical qubit: "
                << this->mapping.getCircQubit(toMoveHwQubit) << '\n';
    } else {
      std::cout << "  not mapped" << '\n';
    }
  }
  nMoves++;
}

qc::Swap qc::NeutralAtomMapper::findBestSwap() {
  // compute necessary movements
  initSwapAndMove(this->frontLayerGate, this->swapCloseByFront,
                  this->moveExactFront);
  initSwapAndMove(this->lookaheadLayerGate, this->swapCloseByLookahead,
                  this->moveExactLookahead);

  // evaluate swaps based on cost function
  auto swaps = getAllPossibleSwaps();
  if (swaps.empty()) {
    return {std::numeric_limits<Qubit>::max(),
            std::numeric_limits<Qubit>::max()};
  }
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
    // return empty swap
    return {std::numeric_limits<Qubit>::max(),
            std::numeric_limits<Qubit>::max()};
  }
  return bestSwap->first;
}

std::set<qc::Swap> qc::NeutralAtomMapper::getAllPossibleSwaps() {
  std::set<Swap> swaps;
  for (const auto& op : this->frontLayerGate) {
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
  auto distanceChangeFront =
      distancePerLayer(swap, this->swapCloseByFront, this->moveExactFront) /
      static_cast<fp>(this->frontLayerGate.size());
  fp distanceChangeLookahead = 0;
  if (!this->lookaheadLayerGate.empty()) {
    distanceChangeLookahead = distancePerLayer(swap, this->swapCloseByLookahead,
                                               this->moveExactLookahead) /
                              static_cast<fp>(this->lookaheadLayerGate.size());
  }
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

void NeutralAtomMapper::initSwapAndMove(
    const qc::GateList& layer, std::vector<SwapOrMove>& swapCloseBy,
    std::vector<std::pair<SwapOrMove, fp>>& moveExact) {
  swapCloseBy.clear();
  moveExact.clear();
  for (const auto& gate : layer) {
    auto usedQubits   = gate->get()->getUsedQubits();
    auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
    if (usedQubits.size() == 2) {
      swapCloseBy.emplace_back(*usedHwQubits.begin(), *usedHwQubits.rbegin());
    } else {
      auto bestPos = getBestMultiQubitPosition(gate);
      if (verbose) {
        std::cout << "bestPos: ";
        for (auto qubit : bestPos) {
          std::cout << qubit << " ";
        }
        std::cout << '\n';
      }
      auto exactMoves = getExactMoveToPosition(gate, bestPos);
      moveExact.insert(moveExact.end(), exactMoves.begin(), exactMoves.end());
    }
  }
}

fp NeutralAtomMapper::distancePerLayer(
    const qc::Swap& swap, const std::vector<SwapOrMove>& swapCloseBy,
    const std::vector<std::pair<SwapOrMove, fp>>& moveExact) {
  fp distBefore = 0;
  fp distAfter  = 0;
  fp distChange = 0;
  // bring close only ontil swap distance =0, bring exact to the exact position
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
  for (const auto& moveWithCost : moveExact) {
    auto move        = moveWithCost.first;
    auto origin      = move.first;
    auto destination = move.second;
    auto cost        = moveWithCost.second;
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
    } else {
      continue;
    }
    distChange +=
        (distAfter - distBefore) / cost * parameters.multiQubitGateWeight;
  }

  return distChange;
}

HwQubits NeutralAtomMapper::getBestMultiQubitPosition(
    std::unique_ptr<qc::Operation>* op) {
  // try to find position around gate Qubits recursively
  // if not, search through coupling graph until found according to a
  // priority queue based on the distance to the other qubits
  std::priority_queue<std::pair<fp, HwQubit>,
                      std::vector<std::pair<fp, HwQubit>>,
                      std::greater<std::pair<fp, HwQubit>>>
       qubitQueue;
  auto gateQubits = op->get()->getUsedQubits();
  // add the gate qubits to the priority queue
  for (const auto& gateQubit : gateQubits) {
    fp totalDist = 0;
    for (const auto& otherGateQubit : gateQubits) {
      if (gateQubit == otherGateQubit) {
        continue;
      }
      totalDist += this->hardwareQubits.getSwapDistance(gateQubit,
                                                        otherGateQubit, false);
    }
    qubitQueue.emplace(totalDist, gateQubit);
  }

  // run through the priority queue until a position is found
  std::set<HwQubit> visitedQubits;
  while (!qubitQueue.empty()) {
    auto qubit = qubitQueue.top().second;
    visitedQubits.insert(qubit);
    qubitQueue.pop();

    auto bestPos = getBestMultiQubitPositionRec(
        gateQubits, {qubit}, this->hardwareQubits.getNearbyQubits(qubit));
    if (!bestPos.empty()) {
      return bestPos;
    }
    // add nearby qubits to the priority queue
    for (const auto& nearbyQubit :
         this->hardwareQubits.getNearbyQubits(qubit)) {
      if (visitedQubits.find(nearbyQubit) != visitedQubits.end()) {
        continue;
      }
      // compute total distance to all other gate qubits
      fp totalDist = 0;
      for (const auto& otherGateQubit : gateQubits) {
        if (nearbyQubit == otherGateQubit) {
          continue;
        }
        totalDist +=
            this->hardwareQubits.getSwapDistance(nearbyQubit, otherGateQubit);
      }
      qubitQueue.emplace(totalDist, nearbyQubit);
    }
  }
  // find gate and move it to the shuttling layer
  if (this->frontLayerGate.find(op) != this->frontLayerGate.end()) {
    this->frontLayerGate.erase(op);
    this->frontLayerShuttling.insert(op);
  }
  // remove from lookahead layer if there
  if (this->lookaheadLayerGate.find(op) != this->lookaheadLayerGate.end()) {
    this->lookaheadLayerGate.erase(op);
    this->lookaheadLayerShuttling.insert(op);
  }
  return {};
}

HwQubits
NeutralAtomMapper::getBestMultiQubitPositionRec(const qc::HwQubits& gateQubits,
                                                qc::Qubits selectedQubits,
                                                qc::Qubits remainingQubits) {
  // if not enough space
  if (selectedQubits.size() + remainingQubits.size() < gateQubits.size()) {
    return {};
  }

  std::vector<std::pair<HwQubit, fp>> summedDistances;
  for (const auto& hwQubit : remainingQubits) {
    fp distance = 0;
    for (const auto& gateHwQubit : gateQubits) {
      distance +=
          this->hardwareQubits.getSwapDistance(hwQubit, gateHwQubit, false);
    }
    summedDistances.emplace_back(hwQubit, distance);
  }
  // select next qubit as the one with minimal distance
  auto nextQubitDist =
      std::min_element(summedDistances.begin(), summedDistances.end(),
                       [](const auto& qubit1, const auto& qubit2) {
                         return qubit1.second < qubit2.second;
                       });
  auto nextQubit = nextQubitDist->first;
  selectedQubits.insert(nextQubit);
  // done if already sufficient number of qubits
  if (selectedQubits.size() == gateQubits.size()) {
    return selectedQubits;
  }
  auto nearbyNextQubit = this->hardwareQubits.getNearbyQubits(nextQubit);
  // compute remaining qubits as the intersection with current
  Qubits newRemainingQubits;
  std::set_intersection(
      remainingQubits.begin(), remainingQubits.end(), nearbyNextQubit.begin(),
      nearbyNextQubit.end(),
      std::inserter(newRemainingQubits, newRemainingQubits.begin()));
  return getBestMultiQubitPositionRec(gateQubits, selectedQubits,
                                      newRemainingQubits);
}

std::vector<std::pair<SwapOrMove, fp>>
NeutralAtomMapper::getExactMoveToPosition(std::unique_ptr<qc::Operation>* op,
                                          HwQubits position) {
  if (position.empty()) {
    return {};
  }
  auto gateQubits = op->get()->getUsedQubits();
  std::vector<std::pair<SwapOrMove, fp>> exactMoves;
  while (!position.empty()) {
    std::vector<std::tuple<HwQubit, std::set<HwQubit>, fp>> minimalDistances;
    std::set<HwQubit> minimalDistancePosQubit;
    for (const auto& gateQubit : gateQubits) {
      fp minimalDistance = std::numeric_limits<fp>::infinity();
      for (const auto& posQubit : position) {
        auto distance =
            this->hardwareQubits.getSwapDistance(gateQubit, posQubit, false);
        if (distance < minimalDistance) {
          minimalDistance = distance;
          minimalDistancePosQubit.clear();
          minimalDistancePosQubit.insert(posQubit);
        } else if (distance == minimalDistance) {
          minimalDistancePosQubit.insert(posQubit);
        }
      }
      if (minimalDistance == std::numeric_limits<fp>::infinity()) {
        // not possible to move to position
        // move gate to shuttling layer
        if (this->frontLayerGate.find(op) != this->frontLayerGate.end()) {
          this->frontLayerGate.erase(op);
          this->frontLayerShuttling.insert(op);
        }
        // remove from lookahead layer if there
        if (this->lookaheadLayerGate.find(op) !=
            this->lookaheadLayerGate.end()) {
          this->lookaheadLayerGate.erase(op);
          this->lookaheadLayerShuttling.insert(op);
        }
        return {};
      }
      minimalDistances.emplace_back(gateQubit, minimalDistancePosQubit,
                                    minimalDistance);
    }
    // find gate qubit with maximal minimal distance to assign first to a
    // position
    auto assignFirst =
        std::max_element(minimalDistances.begin(), minimalDistances.end(),
                         [](const auto& qubit1, const auto& qubit2) {
                           return std::get<2>(qubit1) < std::get<2>(qubit2);
                         });

    auto assignedGateQubit = std::get<0>(*assignFirst);
    auto assignedPosQubits = std::get<1>(*assignFirst);
    // for multiple equal good positions, choose the one that
    // is not assigned to one of the other ones
    HwQubit assignedPosQubit = *assignedPosQubits.begin();
    if (assignedPosQubits.size() > 1) {
      for (const auto& posQubit : assignedPosQubits) {
        // as all places within the position can reach each other, it is
        // sufficient to check for a single unoccupied position
        // check if posQubit is assigned at its current position
        if (std::none_of(minimalDistances.begin(), minimalDistances.end(),
                         [&posQubit](const auto& qubit) {
                           return std::get<0>(qubit) == posQubit &&
                                  *(std::get<1>(qubit).begin()) == posQubit;
                         })) {
          assignedPosQubit = posQubit;
          break;
        }
      }
    }

    // assign gateQubit to position by removing both from gateQubits and
    // position
    gateQubits.erase(gateQubits.find(assignedGateQubit));
    position.erase(position.find(assignedPosQubit));
    // and add to exactMove if not swap with one of the other qubits
    // only problem if their exact swap distance is 1
    if (std::none_of(gateQubits.begin(), gateQubits.end(),
                     [&assignedGateQubit, this](const auto& qubit) {
                       return assignedGateQubit == qubit &&
                              this->hardwareQubits.getSwapDistance(
                                  assignedGateQubit, qubit, false) == 1;
                     }) &&
        assignedGateQubit != assignedPosQubit) {
      exactMoves.emplace_back(
          std::make_pair(assignedGateQubit, assignedPosQubit), 0);
    }
  }

  // compute total distance of all moves
  fp totalDistance = 0;
  for (const auto& [move, cost] : exactMoves) {
    auto [q1, q2] = move;
    totalDistance += this->hardwareQubits.getSwapDistance(q1, q2, false);
  }
  // add cost to the moves -> move first qubit corresponding to almost finished
  // positions
  for (auto& move : exactMoves) {
    // TODO maybe add distance here also to cost -> qubits that are far of
    // are more likely to be moved first
    move.second =
        totalDistance; //- this->hardwareQubits.getSwapDistance(move.first,
                       // move.second, false);
  }

  return exactMoves;
}

AtomMove NeutralAtomMapper::findBestAtomMove() {
  auto moveCombs = getAllPossibleMoveCombinations();

  std::vector<std::pair<AtomMove, fp>> moveCosts;
  moveCosts.reserve(moveCombs.size());
  for (const auto& moveComb : moveCombs) {
    moveCosts.emplace_back(moveComb.getFirstMove(), moveCostComb(moveComb));
  }

  // get move of minimal cost
  auto bestMove = std::min_element(moveCosts.begin(), moveCosts.end(),
                                   [](const auto& move1, const auto& move2) {
                                     return move1.second < move2.second;
                                   });
  return bestMove->first;
}

fp NeutralAtomMapper::moveCostComb(const qc::MoveComb& moveComb) {
  if (!std::isnan(moveComb.cost)) {
    return moveComb.cost;
  }
  fp costComb = 0;
  for (const auto& move : moveComb.moves) {
    costComb += moveCost(move);
  }
  return costComb;
}

fp NeutralAtomMapper::moveCost(const AtomMove& move) {
  fp   cost      = 0;
  auto frontCost = moveDistancePerLayer(move, this->frontLayerShuttling) /
                   static_cast<fp>(this->frontLayerShuttling.size());
  cost += frontCost;
  if (!lookaheadLayerShuttling.empty()) {
    auto lookaheadCost =
        moveDistancePerLayer(move, this->lookaheadLayerShuttling) /
        static_cast<fp>(this->lookaheadLayerShuttling.size());
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
          auto dist    = this->hardwareQubits.getSwapDistance(toMoveHwQubit,
                                                              hwQubit, true);
          distanceBefore += dist;
          if (dist == 0) {
            // add bonus
            distanceBefore -= parameters.shuttlingMakeExecutableBonus *
                              parameters.multiQubitGateWeightShuttling;
          }
        }
        fp distanceAfter = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          auto dist =
              this->hardwareQubits.getSwapDistanceMove(move.second, hwQubit);
          distanceAfter += dist;
          if (dist == 0) {
            // add bonus
            distanceAfter -= parameters.shuttlingMakeExecutableBonus *
                             parameters.multiQubitGateWeightShuttling;
          }
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

MoveCombs NeutralAtomMapper::getAllPossibleMoveCombinations() {
  MoveCombs moves;
  for (const auto& op : this->frontLayerShuttling) {
    auto usedQubits   = op->get()->getUsedQubits();
    auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
    if (usedHwQubits.size() == 2) {
      auto q1      = *usedHwQubits.begin();
      auto q2      = *usedHwQubits.rbegin();
      auto moves12 = getNearbyMoveCombinations(q1, q2);
      auto moves21 = getNearbyMoveCombinations(q2, q1);
      moves.addMoveCombs(moves12);
      moves.addMoveCombs(moves21);
      // TODO test if one should always include move away in both directions
      if (moves.empty()) {
        auto moves1 = getMoveAwayCombinationsNearby(q1, q2);
        auto moves2 = getMoveAwayCombinationsNearby(q2, q1);
        moves.addMoveCombs(moves1);
        moves.addMoveCombs(moves2);
      }
      auto freePos      = getClosestFreePosition(usedQubits);
      auto freePosMoves = getMoveCombinationsToPosition(usedHwQubits, freePos);
      moves.addMoveCombs(freePosMoves);
    } else {
      std::vector<HwQubits> nearbyPositions;
      for (const auto& qubit : usedHwQubits) {
        auto nearbyQubits = this->hardwareQubits.getNearbyQubits(qubit);
        nearbyPositions.push_back(nearbyQubits);
      }
      auto multiQubitMoves =
          getMoveCombinationsToPosition(usedHwQubits, nearbyPositions);
      moves.addMoveCombs(multiQubitMoves);
      auto freePos      = getClosestFreePosition(usedQubits);
      auto freePosMoves = getMoveCombinationsToPosition(usedHwQubits, freePos);
      moves.addMoveCombs(freePosMoves);
    }
  }
  return moves;
}

std::vector<std::set<CoordIndex>>
NeutralAtomMapper::getClosestFreePosition(const qc::Qubits& qubits) {
  // returns a set of coordinates that are free and minimal
  // distance to all other qubits
  auto numFreePos = qubits.size();

  // priority queue ordered by distance of the CoordIndex to
  // all qubits summed up
  std::priority_queue<CoordIndex, std::vector<CoordIndex>,
                      std::function<bool(const CoordIndex&, const CoordIndex&)>>
      q([this, &qubits](const CoordIndex& coord1, const CoordIndex& coord2) {
        fp dist1 = 0;
        fp dist2 = 0;
        for (const auto& qubit : qubits) {
          auto hwQubit = this->mapping.getHwQubit(qubit);
          dist1 += this->arch.getEuclidianDistance(coord1, hwQubit);
          dist2 += this->arch.getEuclidianDistance(coord2, hwQubit);
        }
        return dist1 > dist2;
      });

  for (const auto& qubit : qubits) {
    auto coord = this->hardwareQubits.getCoordIndex(qubit);
    q.push(coord);
  }

  std::vector<CoordIndex> visited;
  while (!q.empty()) {
    auto coord = q.top();
    q.pop();
    visited.push_back(coord);

    auto nearbyFreeCoords =
        this->hardwareQubits.getNearbyFreeCoordinatesByCoord(coord);
    if (nearbyFreeCoords.size() < numFreePos) {
      for (const auto& freeCoord : nearbyFreeCoords) {
        if (std::find(visited.begin(), visited.end(), freeCoord) ==
            visited.end()) {
          q.push(freeCoord);
        }
      }
    } else {
      return {nearbyFreeCoords};
    }
  }
  return {};
}

MoveCombs NeutralAtomMapper::getMoveCombinationsToPosition(
    qc::HwQubits& gateQubits, std::vector<std::set<CoordIndex>>& positions) {
  // compute for each qubit the best position around it based on the cost of the
  // single move choose best one
  if (positions.empty()) {
    return {};
  }
  MoveCombs               moveCombinations;
  std::vector<CoordIndex> gateQubitCoords;
  for (const auto& gateQubit : gateQubits) {
    gateQubitCoords.push_back(this->hardwareQubits.getCoordIndex(gateQubit));
  }

  for (const auto& position : positions) {
    // compute cost for each candidate and each gateQubit
    for (const auto& candidate : position) {
      if (std::find(gateQubitCoords.begin(), gateQubitCoords.end(),
                    candidate) != gateQubitCoords.end()) {
        continue;
      }
      for (const auto& gateQubitCoord : gateQubitCoords) {
        if (gateQubitCoord == candidate ||
            position.find(gateQubitCoord) != position.end()) {
          continue;
        }
        if (this->hardwareQubits.isMapped(candidate)) {
          // add move away move
          auto moveAway = getMoveAwayCombinations(gateQubitCoord, candidate);
          moveCombinations.addMoveCombs(moveAway);
        } else { // candidate coord is free
          auto move = AtomMove{gateQubitCoord, candidate};
          moveCombinations.addMoveComb(MoveComb({move}));
        }
      }
    }
  }

  return moveCombinations;
}

MoveCombs NeutralAtomMapper::getNearbyMoveCombinations(HwQubit start,
                                                       HwQubit target) {
  // Finds free coords near target and moves start to them
  auto      startCoord = this->hardwareQubits.getCoordIndex(start);
  MoveCombs moveCombinations;
  auto nearbyFreeCoords = this->hardwareQubits.getNearbyFreeCoordinates(target);
  for (const auto& coord : nearbyFreeCoords) {
    const AtomMove move = {startCoord, coord};
    moveCombinations.addMoveComb(MoveComb({move}));
  }
  return moveCombinations;
}

MoveCombs NeutralAtomMapper::getMoveAwayCombinationsNearby(qc::HwQubit start,
                                                           qc::HwQubit target) {
  // Finds all possibilities to move a nearby atom away from target
  // and then move start to the now free coord
  MoveCombs moveCombinations;

  auto startCoord         = this->hardwareQubits.getCoordIndex(start);
  auto nearbyTargetQubits = this->hardwareQubits.getNearbyQubits(target);
  for (const auto& targetQubit : nearbyTargetQubits) {
    auto targetCoord = this->hardwareQubits.getCoordIndex(targetQubit);
    auto moveAwayCombinations =
        getMoveAwayCombinations(startCoord, targetCoord);
    moveCombinations.addMoveCombs(moveAwayCombinations);
  }
  return moveCombinations;
}

MoveCombs NeutralAtomMapper::getMoveAwayCombinations(CoordIndex startCoord,
                                                     CoordIndex targetCoord) {
  MoveCombs  moveCombinations;
  auto const originalVector    = this->arch.getVector(startCoord, targetCoord);
  auto const originalDirection = originalVector.direction;
  auto       moveAwayTargets =
      this->hardwareQubits.findClosestFreeCoord(targetCoord, originalDirection);
  for (const auto& moveAwayTarget : moveAwayTargets) {
    const AtomMove move     = {startCoord, targetCoord};
    const AtomMove moveAway = {targetCoord, moveAwayTarget};
    moveCombinations.addMoveComb(MoveComb({moveAway, move}));
  }
  return moveCombinations;
}

std::pair<uint32_t, fp> NeutralAtomMapper::estimateNumSwapGates(
    const std::unique_ptr<qc::Operation>* opPointer) {
  auto* op           = opPointer->get();
  auto  usedQubits   = op->getUsedQubits();
  auto  usedHwQubits = this->mapping.getHwQubits(usedQubits);
  fp    minDistance  = std::numeric_limits<fp>::infinity();
  for (const auto& hwQubit : usedHwQubits) {
    for (const auto& otherHwQubit : usedHwQubits) {
      if (hwQubit == otherHwQubit) {
        continue;
      }
      auto distance =
          this->hardwareQubits.getSwapDistance(hwQubit, otherHwQubit);
      if (distance < minDistance) {
        minDistance = distance;
      }
    }
  }
  const auto minNumSwaps = static_cast<uint32_t>(std::ceil(minDistance));
  const fp   minTime     = minNumSwaps * this->arch.getGateTime(OpType::SWAP);
  return {minNumSwaps, minTime};
}

std::pair<uint32_t, fp> NeutralAtomMapper::estimateNumMove(
    const std::unique_ptr<qc::Operation>* opPointer) {
  auto* op           = opPointer->get();
  auto  usedQubits   = op->getUsedQubits();
  auto  usedHwQubits = this->mapping.getHwQubits(usedQubits);
  auto  usedCoords   = this->hardwareQubits.getCoordIndices(usedHwQubits);
  // estimate the number of moves as:
  // compute distance between qubits
  // 1. for each free coord in the vecinity = 1 move with corresponding distance
  // 2. for each occupied coord in the vecinity = 2 moves with corresponding
  // distance

  uint32_t minMoves = std::numeric_limits<uint32_t>::max();
  fp       minTime  = std::numeric_limits<fp>::infinity();
  for (const auto& coord : usedCoords) {
    fp       totalTime  = 0;
    uint32_t totalMoves = 0;
    auto     nearbyFreeCoords =
        this->hardwareQubits.getNearbyFreeCoordinatesByCoord(coord);
    auto neabyOccupiedCoords =
        this->hardwareQubits.getNearbyOccupiedCoordinatesByCoord(coord);
    auto otherQubitsIt = usedCoords.begin();
    auto nearbyFreeIt  = nearbyFreeCoords.begin();
    auto nearbyOccIt   = neabyOccupiedCoords.begin();
    while (otherQubitsIt != usedCoords.end()) {
      auto otherCoord = *otherQubitsIt;
      if (otherCoord == coord) {
        otherQubitsIt++;
        continue;
      }
      if (nearbyFreeIt != nearbyFreeCoords.end()) {
        totalTime += this->arch.getVectorShuttlingTime(
            this->arch.getVector(otherCoord, *nearbyFreeIt));
        nearbyFreeIt++;
        totalMoves++;
      } else if (nearbyOccIt != neabyOccupiedCoords.end()) {
        totalTime += 2 * this->arch.getVectorShuttlingTime(
                             this->arch.getVector(otherCoord, *nearbyOccIt));
        nearbyOccIt++;
        totalMoves += 2;
      } else {
        throw std::runtime_error(
            "Not enough free or occupied coords to execute a multi-qubit gate. "
            "-> should not happen");
      }

      otherQubitsIt++;
    }

    if (totalTime < minTime) {
      minTime  = totalTime;
      minMoves = totalMoves;
    }
  }

  return {minMoves, minTime};
}

bool NeutralAtomMapper::swapGateBetter(
    const std::unique_ptr<qc::Operation>* opPointer) {
  auto [minNumSwaps, minTimeSwaps] = estimateNumSwapGates(opPointer);
  auto [minMoves, minTimeMoves]    = estimateNumMove(opPointer);
  auto fidSwaps =
      std::exp(-minTimeSwaps * this->arch.getNqubits() /
               this->arch.getDecoherenceTime()) *
      std::pow(this->arch.getGateAverageFidelity(OpType::SWAP), minNumSwaps);
  auto fidMoves =
      std::exp(-minTimeMoves * this->arch.getNqubits() /
               this->arch.getDecoherenceTime()) *
      std::pow(this->arch.getGateAverageFidelity(OpType::Move), minMoves);

  return fidSwaps * parameters.gateWeight >
         fidMoves * parameters.shuttlingWeight;
}

} // namespace qc
