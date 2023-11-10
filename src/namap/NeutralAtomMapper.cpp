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
      std::vector<const Operation*> gatesToExecute;
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
      GateList gatesToExecute;
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
  return mappedQc;
}

QuantumComputation NeutralAtomMapper::mapAod(qc::QuantumComputation& qc) {
  // decompose SWAP gates
  CircuitOptimizer::decomposeSWAP(qc, false);
  CircuitOptimizer::replaceMCXWithMCZ(qc);
  CircuitOptimizer::singleQubitGateFusion(qc);
  CircuitOptimizer::flattenOperations(qc);
  // decompose AOD moves
  AodScheduler scheduler(this->arch);
  auto         aodQc = scheduler.schedule(qc);
  std::cout << "nMoveGroups: " << scheduler.getNMoveGroups() << '\n';
  return aodQc;
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

void qc::NeutralAtomMapper::updateFrontLayerByGate(GateList& gatesToExecute) {
  for (const auto& gate : gatesToExecute) {
    // add other gate qubits to be checked
    for (auto qubit : gate->getUsedQubits()) {
      this->frontQubitsToUpdate.insert(qubit);
      this->lookaheadOffsets[qubit]--;
      this->frontLayerIterators[qubit]++;
    }
    // remove from current FrontLayer
    mapGate(gate);
    // update DAG iterators
    // remove from FrontLayer
    if (std::find(this->frontLayerGate.begin(), this->frontLayerGate.end(),
                  gate) != this->frontLayerGate.end()) {
      this->frontLayerGate.erase(std::find(this->frontLayerGate.begin(),
                                           this->frontLayerGate.end(), gate));
    }
    if (std::find(this->frontLayerShuttling.begin(),
                  this->frontLayerShuttling.end(),
                  gate) != this->frontLayerShuttling.end()) {
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
  for (const auto& qubit : this->frontQubitsToUpdate) {
    auto tempIter = this->frontLayerIterators[qubit];
    while (tempIter < this->dag[qubit].end()) {
      auto* op = (*tempIter)->get();
      if (op->getUsedQubits().size() == 1) {
        mapGate(op);
        this->frontLayerIterators[qubit]++;
        tempIter++;
      } else {
        // continue if following gates commute
        bool commutes = true;
        while (commutes && tempIter < this->dag[qubit].end()) {
          auto* nextOp = (*tempIter)->get();
          commutes =
              commutesWithAtQubit(this->frontLayerGate, nextOp, qubit) &&
              commutesWithAtQubit(this->frontLayerShuttling, nextOp, qubit) &&
              commutesWithAtQubit(this->frontCandidates[qubit], nextOp, qubit);
          if (commutes) {
            if (nextOp->getUsedQubits().size() == 1) {
              mapGate(nextOp);
            } else { // not executable but commutes
              if (std::find(this->frontLayerGate.begin(),
                            this->frontLayerGate.end(),
                            nextOp) == this->frontLayerGate.end() &&
                  std::find(this->frontLayerShuttling.begin(),
                            this->frontLayerShuttling.end(),
                            nextOp) == this->frontLayerShuttling.end()) {
                this->frontCandidates[qubit].push_back(nextOp);
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
  for (const auto& qubit : tempQubitsToUpdate) {
    std::vector<const Operation*> toRemove;
    for (auto* opPointer : this->frontCandidates[qubit]) {
      bool inFrontLayer = true;
      if (opPointer->getType() == qc::OpType::I) {
        toRemove.push_back(opPointer);
        continue;
      }
      for (const auto& opQubit : opPointer->getUsedQubits()) {
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
          for (const auto& opQubit : opPointer->getUsedQubits()) {
            this->frontLayerIterators[opQubit]++;
            this->frontQubitsToUpdate.insert(opQubit);
          }
          // remove from lookahead layer if there
          auto inLookaheadGate =
              std::find(this->lookaheadLayerGate.begin(),
                        this->lookaheadLayerGate.end(), opPointer);
          if (inLookaheadGate != this->lookaheadLayerGate.end()) {
            this->lookaheadLayerGate.erase(inLookaheadGate);
          }
          auto inLookaheadShuttling =
              std::find(this->lookaheadLayerShuttling.begin(),
                        this->lookaheadLayerShuttling.end(), opPointer);
          if (inLookaheadShuttling != this->lookaheadLayerShuttling.end()) {
            this->lookaheadLayerShuttling.erase(inLookaheadShuttling);
          }
        } else { // not executable
          addToFrontLayer(opPointer);
          // move lookahead
          for (const auto& opQubit : opPointer->getUsedQubits()) {
            this->lookaheadOffsets[opQubit]++;
          }
        }
        // remove from candidacy of other qubits
        for (const auto& opQubit : opPointer->getUsedQubits()) {
          this->lookaheadQubitsToUpdate.insert(opQubit);
          if (qubit == opQubit) {
            continue;
          }
          this->frontCandidates[opQubit].erase(
              std::find(this->frontCandidates[opQubit].begin(),
                        this->frontCandidates[opQubit].end(), opPointer));
        }

        // save to remove from candidacy of this qubit
        toRemove.push_back(opPointer);
      }
    }
    // remove from candidacy of this qubit
    // has to be done now to not change iterating list
    for (const auto* opPointer : toRemove) {
      this->frontCandidates[qubit].erase(
          std::find(this->frontCandidates[qubit].begin(),
                    this->frontCandidates[qubit].end(), opPointer));
    }
  }
  // if gates have been executed directly refresh front layer
  if (!this->frontQubitsToUpdate.empty()) {
    updateFrontLayerByQubit();
  }
}

bool qc::NeutralAtomMapper::commutesWithAtQubit(const qc::GateList& layer,
                                                const Operation*    opPointer,
                                                const Qubit&        qubit) {
  return std::all_of(layer.begin(), layer.end(),
                     [&opPointer, &qubit](const auto& frontOpPointer) {
                       return commuteAtQubit(opPointer, frontOpPointer, qubit);
                     });
}

bool qc::NeutralAtomMapper::commuteAtQubit(const Operation* op1,
                                           const Operation* op2,
                                           const qc::Qubit& qubit) {
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
  // commute if both are controlled at qubit or const Operation* on qubit is
  // same check controlles
  if (op1->getControls().find(qubit) != op1->getControls().end() &&
      op2->getControls().find(qubit) != op2->getControls().end()) {
    return true;
  }
  // control and Z also commute
  if ((op1->getControls().find(qubit) != op1->getControls().end() &&
       op2->getType() == qc::OpType::Z) ||
      (op2->getControls().find(qubit) != op2->getControls().end() &&
       op1->getType() == qc::OpType::Z)) {
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
    std::vector<const Operation*> toRemove;
    for (auto* opPointer : this->lookaheadCandidates[qubit]) {
      bool inLookahead = true;
      for (const auto& opQubit : opPointer->getUsedQubits()) {
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
        for (const auto& opQubit : opPointer->getUsedQubits()) {
          if (qubit == opQubit) {
            continue;
          }
          this->lookaheadCandidates[opQubit].erase(
              std::find(this->lookaheadCandidates[opQubit].begin(),
                        this->lookaheadCandidates[opQubit].end(), opPointer));
        }
        // save to remove from candidacy of this qubit
        toRemove.push_back(opPointer);
      }
    }
    // remove from candidacy of this qubit
    for (const auto* opPointer : toRemove) {
      this->lookaheadCandidates[qubit].erase(
          std::find(this->lookaheadCandidates[qubit].begin(),
                    this->lookaheadCandidates[qubit].end(), opPointer));
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
    while (lookaheadIter < this->dag[qubit].end() &&
           foundMultiQubitGate < this->lookaheadDepth) {
      auto* opPointer = (*lookaheadIter)->get();
      if (opPointer->getType() == qc::OpType::I) {
        lookaheadIter++;
        continue;
      }
      if (opPointer->getUsedQubits().size() > 1) {
        foundMultiQubitGate++;
        if (std::find(this->lookaheadCandidates[qubit].begin(),
                      this->lookaheadCandidates[qubit].end(),
                      opPointer) == this->lookaheadCandidates[qubit].end()) {
          this->lookaheadCandidates[qubit].push_back(opPointer);
        }
        lookaheadIter++;
        // continue if following gates commute
        // only check last gate (only approximate)
        bool commutes = true;
        while (commutes && lookaheadIter < this->dag[qubit].end()) {
          auto* nextOpPointer = (*lookaheadIter)->get();
          commutes            = commuteAtQubit(opPointer, nextOpPointer, qubit);
          if (commutes) {
            lookaheadIter++;
            if (nextOpPointer->getUsedQubits().size() > 1) {
              if (std::find(this->lookaheadCandidates[qubit].begin(),
                            this->lookaheadCandidates[qubit].end(),
                            nextOpPointer) ==
                  this->lookaheadCandidates[qubit].end()) {
                this->lookaheadCandidates[qubit].push_back(nextOpPointer);
              }
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

void qc::NeutralAtomMapper::mapGate(const Operation* op) {
  if (op->getType() == qc::OpType::I) {
    return;
  }
  if (this->verbose) {
    std::cout << "mapped " << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << "\n";
  }
  // convert circuit qubits to CoordIndex and append to mappedQc
  auto  opCopyUnique = op->clone();
  auto* opCopy       = opCopyUnique.get();
  this->mapping.mapToHwQubits(opCopy);
  this->hardwareQubits.mapToCoordIdx(opCopy);
  this->mappedQc.emplace_back(opCopy->clone());
}

bool qc::NeutralAtomMapper::isExecutable(const Operation* opPointer) {
  auto usedQubits  = opPointer->getUsedQubits();
  auto nUsedQubits = usedQubits.size();
  if (nUsedQubits == 1) {
    return true;
  }
  std::set<Qubit> usedHwQubits;
  for (auto qubit : usedQubits) {
    usedHwQubits.insert(this->mapping.getHwQubit(qubit));
  }
  return this->hardwareQubits.getTotalDistance(usedHwQubits) == 0;
}

void qc::NeutralAtomMapper::addToFrontLayer(const Operation* opPointer) {
  if (swapGateBetter(opPointer)) {
    this->frontLayerGate.push_back(opPointer);
    // remove from lookahead layer if there
    auto idxLookaheadGate =
        std::find(this->lookaheadLayerGate.begin(),
                  this->lookaheadLayerGate.end(), opPointer);
    if (idxLookaheadGate != this->lookaheadLayerGate.end()) {
      this->lookaheadLayerGate.erase(idxLookaheadGate);
    }
  } else {
    this->frontLayerShuttling.push_back(opPointer);
    // remove from lookahead layer if there
    auto idxLookaheadShuttling =
        std::find(this->lookaheadLayerShuttling.begin(),
                  this->lookaheadLayerShuttling.end(), opPointer);
    if (idxLookaheadShuttling != this->lookaheadLayerShuttling.end()) {
      this->lookaheadLayerShuttling.erase(idxLookaheadShuttling);
    }
  }
  // remove from front candidates not allowed here
  // to prevent modify while iterating
}

void NeutralAtomMapper::addToLookaheadLayer(const Operation* opPointer) {
  if (swapGateBetter(opPointer)) {
    this->lookaheadLayerGate.push_back(opPointer);
  } else {
    this->lookaheadLayerShuttling.push_back(opPointer);
  }
}

void qc::NeutralAtomMapper::printLayers() {
  std::cout << "f,g: ";
  for (const auto* op : this->frontLayerGate) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << "f,s: ";
  for (const auto* op : this->frontLayerShuttling) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << "l,g: ";
  for (const auto* op : this->lookaheadLayerGate) {
    std::cout << op->getName() << " ";
    for (auto qubit : op->getUsedQubits()) {
      std::cout << qubit << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
  std::cout << "l,g: ";
  for (const auto* op : this->lookaheadLayerShuttling) {
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
      executableGates.push_back(opPointer);
    }
  }
  for (auto* opPointer : this->frontLayerShuttling) {
    if (isExecutable(opPointer)) {
      executableGates.push_back(opPointer);
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
    for (const auto& qubit : op->getUsedQubits()) {
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
    auto usedQubits   = gate->getUsedQubits();
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

HwQubits NeutralAtomMapper::getBestMultiQubitPosition(const Operation* op) {
  // try to find position around gate Qubits recursively
  // if not, search through coupling graph until found according to a
  // priority queue based on the distance to the other qubits
  std::priority_queue<std::pair<fp, HwQubit>,
                      std::vector<std::pair<fp, HwQubit>>,
                      std::greater<std::pair<fp, HwQubit>>>
       qubitQueue;
  auto gateQubits = op->getUsedQubits();
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
  auto idxFrontGate =
      std::find(this->frontLayerGate.begin(), this->frontLayerGate.end(), op);
  if (idxFrontGate != this->frontLayerGate.end()) {
    this->frontLayerGate.erase(idxFrontGate);
    this->frontLayerShuttling.push_back(op);
  }
  // remove from lookahead layer if there
  auto idxLookaheadGate = std::find(this->lookaheadLayerGate.begin(),
                                    this->lookaheadLayerGate.end(), op);
  if (idxLookaheadGate != this->lookaheadLayerGate.end()) {
    this->lookaheadLayerGate.erase(idxLookaheadGate);
    this->lookaheadLayerShuttling.push_back(op);
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
NeutralAtomMapper::getExactMoveToPosition(const Operation* op,
                                          HwQubits         position) {
  if (position.empty()) {
    return {};
  }
  auto                                   gateQubits = op->getUsedQubits();
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
        auto idxFrontGate = std::find(this->frontLayerGate.begin(),
                                      this->frontLayerGate.end(), op);
        if (idxFrontGate != this->frontLayerGate.end()) {
          this->frontLayerGate.erase(idxFrontGate);
          this->frontLayerShuttling.push_back(op);
        }
        // remove from lookahead layer if there
        auto idxLookaheadGate = std::find(this->lookaheadLayerGate.begin(),
                                          this->lookaheadLayerGate.end(), op);
        if (idxLookaheadGate != this->lookaheadLayerGate.end()) {
          this->lookaheadLayerGate.erase(idxLookaheadGate);
          this->lookaheadLayerShuttling.push_back(op);
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
  auto moveCombs = getBestPossibleMoveCombinations();

  std::vector<std::pair<MoveComb, fp>> moveCosts;
  moveCosts.reserve(moveCombs.size());
  for (const auto& moveComb : moveCombs) {
    moveCosts.emplace_back(moveComb, moveCostComb(moveComb));
  }

  std::sort(moveCosts.begin(), moveCosts.end(),
            [](const auto& move1, const auto& move2) {
              return move1.second < move2.second;
            });

  // get move of minimal cost
  auto bestMove = std::min_element(moveCosts.begin(), moveCosts.end(),
                                   [](const auto& move1, const auto& move2) {
                                     return move1.second < move2.second;
                                   });
  return bestMove->first.getFirstMove();
}

fp NeutralAtomMapper::moveCostComb(const qc::MoveComb& moveComb) {
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
  if (!this->lastMoves.empty()) {
    auto parallelCost = parameters.shuttlingTimeWeight *
                        parallelMoveCost(move) /
                        static_cast<fp>(this->lastMoves.size()) /
                        static_cast<fp>(this->frontLayerShuttling.size());
    cost += parallelCost;
  }
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
      auto usedQubits = gate->getUsedQubits();
      if (usedQubits.find(toMoveCircuitQubit) != usedQubits.end()) {
        // check distance reduction
        fp distanceBefore   = 0;
        fp executableBefore = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          auto dist    = this->hardwareQubits.getSwapDistance(toMoveHwQubit,
                                                              hwQubit, true);
          //                    distanceBefore += dist;
          if (dist == 0) {
            executableBefore++;
          }
        }
        fp distanceAfter   = 0;
        fp executableAfter = 0;
        for (const auto& qubit : usedQubits) {
          if (qubit == toMoveCircuitQubit) {
            continue;
          }
          auto hwQubit = this->mapping.getHwQubit(qubit);
          auto dist =
              this->hardwareQubits.getSwapDistanceMove(move.second, hwQubit);
          //                    distanceAfter += dist;
          if (dist == 0) {
            executableAfter++;
          }
        }
        distChange += distanceAfter - distanceBefore;
        const auto execBonus = (executableAfter - executableBefore) *
                               parameters.shuttlingMakeExecutableBonus;
        distChange += execBonus;
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
  // check if in same row/column like last moves
  // then can may be loaded in parallel
  auto moveCoordInit = this->arch.getCoordinate(move.first);
  auto moveCoordEnd  = this->arch.getCoordinate(move.second);
  parallelCost += arch.getShuttlingTime(OpType::AodActivate) +
                  arch.getShuttlingTime(OpType::AodDeactivate);
  for (const auto& lastMove : this->lastMoves) {
    auto lastMoveCoordInit = this->arch.getCoordinate(lastMove.first);
    auto lastMoveCoordEnd  = this->arch.getCoordinate(lastMove.second);
    if (moveCoordInit.getX() == lastMoveCoordInit.getX() ||
        moveCoordInit.getY() == lastMoveCoordInit.getY()) {
      parallelCost -= arch.getShuttlingTime(OpType::AodActivate);
    }
    if (moveCoordEnd.getX() == lastMoveCoordEnd.getX() ||
        moveCoordEnd.getY() == lastMoveCoordEnd.getY()) {
      parallelCost -= arch.getShuttlingTime(OpType::AodDeactivate);
    }
  }
  // check if move can use AOD atom from last moves
  //  if (std::find(lastEndingCoords.begin(), lastEndingCoords.end(),
  //  move.first) ==
  //      lastEndingCoords.end()) {
  //    parallelCost += arch.getShuttlingTime(OpType::AodActivate) +
  //                    arch.getShuttlingTime(OpType::AodDeactivate);
  //  }
  return parallelCost;
}

MultiQubitMovePos
NeutralAtomMapper::getMovePositionRec(MultiQubitMovePos   currentPos,
                                      const CoordIndices& gateCoords,
                                      const size_t&       maxNMoves) {
  if (currentPos.coords.size() == gateCoords.size()) {
    return currentPos;
  }
  if (currentPos.nMoves > maxNMoves) {
    return {};
  }

  auto nearbyCoords = this->arch.getNearbyCoordinates(currentPos.coords.back());
  // filter out coords that have a SWAP distance unequal to 0 to any of the
  // current qubits. Also sort out coords that are already in the vector
  std::vector<CoordIndex> filteredNearbyCoords;
  for (const auto& coord : nearbyCoords) {
    bool valid = true;
    for (const auto& qubit : currentPos.coords) {
      if (this->arch.getSwapDistance(qubit, coord) != 0 || coord == qubit) {
        valid = false;
        break;
      }
    }
    if (valid) {
      filteredNearbyCoords.push_back(coord);
    }
  }

  // differentiate between free and occupied coords
  CoordIndices freeNearbyCoords;
  CoordIndices occupiedNearbyCoords;
  CoordIndices occupiedGateCoords;
  for (const auto& coord : filteredNearbyCoords) {
    if (this->hardwareQubits.isMapped(coord)) {
      if (std::find(gateCoords.begin(), gateCoords.end(), coord) !=
          gateCoords.end()) {
        occupiedGateCoords.push_back(coord);
      } else {
        occupiedNearbyCoords.push_back(coord);
      }
    } else {
      freeNearbyCoords.push_back(coord);
    }
  }

  // compute minimal possible moves
  size_t       minPossibleMoves = currentPos.nMoves;
  size_t const nMissingQubits   = gateCoords.size() - currentPos.coords.size();
  auto         itGate           = occupiedGateCoords.begin();
  auto         itFree           = freeNearbyCoords.begin();
  auto         itOcc            = occupiedNearbyCoords.begin();
  for (size_t i = 0; i < nMissingQubits; ++i) {
    if (itGate != occupiedGateCoords.end()) {
      ++itGate;
    } else if (itFree != freeNearbyCoords.end()) {
      ++itFree;
      minPossibleMoves += 1;
    } else if (itOcc != occupiedNearbyCoords.end()) {
      ++itOcc;
      minPossibleMoves += 2;
    }
  }
  if (minPossibleMoves > maxNMoves) {
    return {};
  }

  for (const auto& gateCoord : occupiedGateCoords) {
    MultiQubitMovePos nextPos = MultiQubitMovePos(currentPos);
    nextPos.coords.push_back(gateCoord);
    auto bestPos = getMovePositionRec(nextPos, gateCoords, maxNMoves);
    if (bestPos.coords.size() == gateCoords.size()) {
      return bestPos;
    }
  }

  for (const auto& freeCoord : freeNearbyCoords) {
    MultiQubitMovePos nextPos = MultiQubitMovePos(currentPos);
    nextPos.coords.push_back(freeCoord);
    nextPos.nMoves += 1;
    auto bestPos = getMovePositionRec(nextPos, gateCoords, maxNMoves);
    if (bestPos.coords.size() == gateCoords.size()) {
      return bestPos;
    }
  }

  for (const auto& occCoord : occupiedNearbyCoords) {
    MultiQubitMovePos nextPos = MultiQubitMovePos(currentPos);
    nextPos.coords.push_back(occCoord);
    nextPos.nMoves += 2;
    auto bestPos = getMovePositionRec(nextPos, gateCoords, maxNMoves);
    if (bestPos.coords.size() == gateCoords.size()) {
      return bestPos;
    }
  }

  // if no position found, return empty
  return {};
}

MoveCombs NeutralAtomMapper::getBestPossibleMoveCombinations() {
  MoveCombs allMoves;
  for (const auto& op : this->frontLayerShuttling) {
    auto usedQubits    = op->getUsedQubits();
    auto usedHwQubits  = this->mapping.getHwQubits(usedQubits);
    auto usedCoordsSet = this->hardwareQubits.getCoordIndices(usedHwQubits);
    auto usedCoords =
        std::vector<CoordIndex>(usedCoordsSet.begin(), usedCoordsSet.end());
    auto bestPos = getBestMovePos(usedCoords);
    auto moves   = getMoveCombinationsToPosition(usedHwQubits, bestPos.coords);
    allMoves.addMoveCombs(moves);
  }
  allMoves.removeLongerMoveCombs();
  return allMoves;
}

MultiQubitMovePos
NeutralAtomMapper::getBestMovePos(const CoordIndices& gateCoords) {
  size_t const maxMoves   = gateCoords.size() * 2;
  size_t const minMoves   = gateCoords.size();
  size_t       nMovesGate = maxMoves;
  // do a breadth first search for the best position
  // start with the used coords
  std::queue<CoordIndex> q;
  for (const auto& coord : gateCoords) {
    q.push(coord);
  }
  std::vector<CoordIndex> visited;

  auto finalBestPos = MultiQubitMovePos();
  while (!q.empty()) {
    auto coord = q.front();
    q.pop();
    if (std::find(visited.begin(), visited.end(), coord) != visited.end()) {
      continue;
    }
    visited.push_back(coord);
    MultiQubitMovePos currentPos;
    currentPos.coords.push_back(coord);
    if (this->hardwareQubits.isMapped(coord)) {
      if (std::find(gateCoords.begin(), gateCoords.end(), coord) !=
          gateCoords.end()) {
        currentPos.nMoves = 0;
      } else {
        currentPos.nMoves = 2;
      }
    } else {
      currentPos.nMoves = 1;
    }
    auto bestPos = getMovePositionRec(currentPos, gateCoords, nMovesGate);
    finalBestPos = bestPos;
    if (!bestPos.coords.empty() && bestPos.nMoves <= minMoves) {
      return bestPos;
    }

    // min not yet reached, check nearby
    if (!bestPos.coords.empty()) {
      nMovesGate = std::min(nMovesGate, bestPos.nMoves);
    }
    for (const auto& nearbyCoord : this->arch.getNearbyCoordinates(coord)) {
      if (std::find(visited.begin(), visited.end(), nearbyCoord) ==
          visited.end()) {
        q.push(nearbyCoord);
      }
    }
  }
  throw QFRException(
      "No move position found (check if enough free coords are available)");
  return finalBestPos;
}

MoveCombs NeutralAtomMapper::getMoveCombinationsToPosition(
    qc::HwQubits& gateQubits, std::vector<CoordIndex>& position) {
  if (position.empty()) {
    throw QFRException("No position given");
  }
  // compute for each qubit the best position around it based on the cost of
  // the single move choose best one
  MoveCombs            moveCombinations;
  std::set<CoordIndex> gateQubitCoords;
  for (const auto& gateQubit : gateQubits) {
    gateQubitCoords.insert(this->hardwareQubits.getCoordIndex(gateQubit));
  }

  auto     remainingCoords = position;
  MoveComb moveComb;
  // compute cost for each candidate and each gateQubit
  auto remainingGateCoords = gateQubitCoords;
  // pre-filter away all gateQubitCoords which are already in the position
  for (auto it = remainingGateCoords.begin();
       it != remainingGateCoords.end();) {
    if (std::find(remainingCoords.begin(), remainingCoords.end(), *it) !=
        remainingCoords.end()) {
      remainingCoords.erase(
          std::find(remainingCoords.begin(), remainingCoords.end(), *it));
      it = remainingGateCoords.erase(it);
    } else {
      ++it;
    }
  }

  while (!remainingGateCoords.empty()) {
    auto currentGateQubit = *remainingGateCoords.begin();
    // compute costs and find best coord
    std::vector<std::pair<CoordIndex, fp>> costs;
    for (const auto& remainingCoord : remainingCoords) {
      auto cost = moveCost({currentGateQubit, remainingCoord});
      costs.emplace_back(remainingCoord, cost);
    }
    // find minimal cost
    auto bestCost  = std::min_element(costs.begin(), costs.end(),
                                      [](const auto& cost1, const auto& cost2) {
                                       return cost1.second < cost2.second;
                                     });
    auto bestCoord = bestCost->first;
    if (this->hardwareQubits.isMapped(bestCoord)) {
      auto moveAwayComb = getMoveAwayCombinations(currentGateQubit, bestCoord);
      for (const auto& moveAway : moveAwayComb) {
        moveComb.append(moveAway);
      }
    } else { // free coord
      moveComb.append(AtomMove{currentGateQubit, bestCoord});
    }
    remainingGateCoords.erase(currentGateQubit);
    remainingCoords.erase(
        std::find(remainingCoords.begin(), remainingCoords.end(), bestCoord));
  }
  return MoveCombs({moveComb});
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

std::pair<uint32_t, fp>
NeutralAtomMapper::estimateNumSwapGates(const Operation* opPointer) {
  auto usedQubits   = opPointer->getUsedQubits();
  auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
  fp   minDistance  = std::numeric_limits<fp>::infinity();
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
  const fp   minTime     = minNumSwaps * this->arch.getGateTime("swap");
  return {minNumSwaps, minTime};
}

std::pair<uint32_t, fp>
NeutralAtomMapper::estimateNumMove(const Operation* opPointer) {
  auto usedQubits   = opPointer->getUsedQubits();
  auto usedHwQubits = this->mapping.getHwQubits(usedQubits);
  auto usedCoords   = this->hardwareQubits.getCoordIndices(usedHwQubits);
  // estimate the number of moves as:
  // compute distance between qubits
  // 1. for each free coord in the vecinity = 1 move with corresponding
  // distance
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
        totalTime += this->arch.getShuttlingTime(OpType::AodActivate) +
                     this->arch.getShuttlingTime(OpType::AodDeactivate);
        nearbyFreeIt++;
        totalMoves++;
      } else if (nearbyOccIt != neabyOccupiedCoords.end()) {
        totalTime += 2 * this->arch.getVectorShuttlingTime(
                             this->arch.getVector(otherCoord, *nearbyOccIt));
        totalTime += 2 * (this->arch.getShuttlingTime(OpType::AodActivate) +
                          this->arch.getShuttlingTime(OpType::AodDeactivate));
        nearbyOccIt++;
        totalMoves += 2;
      } else {
        throw std::runtime_error("No space to "
                                 "execute a multi-qubit gate. "
                                 "Check int radius. Op:" +
                                 opPointer->getName() + " nQubit: " +
                                 std::to_string(usedQubits.size()));
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

bool NeutralAtomMapper::swapGateBetter(const Operation* opPointer) {
  auto [minNumSwaps, minTimeSwaps] = estimateNumSwapGates(opPointer);
  auto [minMoves, minTimeMoves]    = estimateNumMove(opPointer);
  auto fidSwaps =
      std::exp(-minTimeSwaps * this->arch.getNqubits() /
               this->arch.getDecoherenceTime()) *
      std::pow(this->arch.getGateAverageFidelity("swap"), minNumSwaps);
  auto fidMoves =
      std::exp(-minTimeMoves * this->arch.getNqubits() /
               this->arch.getDecoherenceTime()) *
      std::pow(this->arch.getShuttlingAverageFidelity(OpType::AodMove),
               minMoves);

  return fidSwaps * parameters.gateWeight >
         fidMoves * parameters.shuttlingWeight;
}

} // namespace qc
