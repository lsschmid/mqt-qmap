//
// Created by Ludwig Schmid on 05.10.23.
//

#include "namap/NeutralAtomMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "namap/NeutralAtomUtils.hpp"
#include "utils.hpp"

namespace qc {
void qc::NeutralAtomMapper::map(qc::QuantumComputation& qc) {
  this->dag     = qc::CircuitOptimizer::constructDAG(qc);
  this->mapping = Mapping(qc.getNqubits(), InitialMapping::Identity);
  for (auto& i : dag) {
    this->frontLayerIterators.push_back(i.begin());
    this->lookaheadOffsets.push_back(0);
    this->lookaheadCandidates.emplace_back();
    this->frontCandidates.emplace_back();
  }

  std::cout << "test" << std::endl;

  //   precompute all distances

  createFrontLayer();
  auto gatesToExecute = getExecutableGates();
  updateFrontLayerByGate(gatesToExecute);
  //  printLayers();
  updateLookaheadLayerByQubit();
  printLayers();

  // probably need
  // updateFrontLayerByQubit(list with all qubits) to trigger commutation

  // TODO: find SWAP
  // TODO: get list of front gates that can be executed now
  auto i = 0;
  while (!this->frontLayer.empty()) {
    i++;
    std::cout << "iteration " << i << std::endl;
    std::set<std::unique_ptr<Operation>*> gatesToExecute;
    while (gatesToExecute.empty()) {
      auto bestSwap = findBestSwap();
      updateMapping(bestSwap);
      gatesToExecute = getExecutableGates();
    }
    updateFrontLayerByGate(gatesToExecute);
    gatesToExecute = getExecutableGates();
    updateFrontLayerByGate(gatesToExecute);
    //    printLayers();
    updateLookaheadLayerByQubit();
    printLayers();
  }
}

void qc::NeutralAtomMapper::createFrontLayer() {
  // vector indicating if qubit is free
  std::vector<Qubit> allQubits;
  for (uint32_t i = 0; i < this->dag.size(); i++) {
    allQubits.push_back(i);
  }
  this->lookaheadQubitsToUpdate = allQubits;
  updateLookaheadLayerByQubit();
  this->lookaheadQubitsToUpdate.clear();
  this->frontQubitsToUpdate = allQubits;
  updateFrontLayerByQubit();
}

void qc::NeutralAtomMapper::updateFrontLayerByGate(
    std::set<std::unique_ptr<Operation>*>& gatesToExecute) {
  for (auto* gate : gatesToExecute) {
    // add other gate qubits to be checked
    for (auto qubit : (*gate)->getUsedQubits()) {
      if (std::find(this->frontQubitsToUpdate.begin(),
                    this->frontQubitsToUpdate.end(),
                    qubit) == this->frontQubitsToUpdate.end()) {
        this->frontQubitsToUpdate.push_back(qubit);
      }
    }
    // remove from current FrontLayer
    mapGate(gate);
    // update DAG iterators
    for (auto gateQubit : (*gate)->getUsedQubits()) {
      this->lookaheadOffsets[gateQubit]--;
    }
    // remove from FrontLayer
    this->frontLayer.erase(
        std::find(this->frontLayer.begin(), this->frontLayer.end(), gate));
  }
  updateFrontLayerByQubit();
}

void qc::NeutralAtomMapper::updateFrontLayerByQubit() {
  findFrontCandidates();
  updateFrontLayerByCandidates();
  this->frontQubitsToUpdate.clear();
}

void qc::NeutralAtomMapper::findFrontCandidates() {
  for (auto& qubit : this->frontQubitsToUpdate) {
    auto tempIter = this->frontLayerIterators[qubit];
    while (tempIter != this->dag[qubit].end()) {
      auto* opPointer = *tempIter;
      if (opPointer->get()->getUsedQubits().size() == 1) {
        if (std::find(this->executedCommutingGates.rbegin(),
                      this->executedCommutingGates.rend(),
                      opPointer) == this->executedCommutingGates.rend()) {
          mapGate(opPointer);
        }
        tempIter++;
      } else {
        this->frontCandidates[qubit].insert(opPointer);
        tempIter++;

        // continue if following gates commute
        bool commutes = true;
        while (commutes && tempIter != this->dag[qubit].end()) {
          auto* nextOpPointer = *tempIter;
          commutes            = commutesWith(this->frontLayer, nextOpPointer) &&
                     commutesWith(this->frontCandidates[qubit], nextOpPointer);
          if (commutes) {
            if (nextOpPointer->get()->getUsedQubits().size() == 1) {
              mapGate(nextOpPointer);
              this->executedCommutingGates.push_back(nextOpPointer);
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
  std::vector<Qubit> qubitsToCheckNext;
  for (auto& qubit : this->frontQubitsToUpdate) {
    std::vector<std::unique_ptr<qc::Operation>*> toRemove;
    for (auto* opPointer : this->frontCandidates[qubit]) {
      bool  inFrontLayer = true;
      auto* op           = opPointer->get();
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
        }
        addToFrontLayer(opPointer);
        // remove from candidacy of other qubits
        for (const auto& opQubit : op->getUsedQubits()) {
          qubitsToCheckNext.push_back(opQubit);
          this->lookaheadOffsets[opQubit]++;
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
    for (auto* opPointer : toRemove) {
      this->frontCandidates[qubit].erase(opPointer);
    }
  }
  this->lookaheadQubitsToUpdate = qubitsToCheckNext;
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
  this->mappedQc.emplace_back((*op)->clone());
  for (auto qubit : (*op)->getUsedQubits()) {
    this->frontLayerIterators[qubit]++;
    if (qubit == 0) {
      this->debugCounter++;
      std::cout << "debugCounter single commute " << this->debugCounter
                << std::endl;
    }
  }

  std::cout << "mapped " << (*op)->getName();
  for (auto qubit : (*op)->getUsedQubits()) {
    std::cout << " " << qubit;
  }
  std::cout << std::endl;
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
  if (nUsedQubits == 2) {
    return this->hardwareQubits.getTotalDistance(usedHwQubits) == 0;
  }
  return false;
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
  this->mapping.swap(swap);
  this->mappedQc.swap(swap.first, swap.second);
  std::cout << "swapped " << swap.first << " " << swap.second << std::endl;
}

qc::Swap qc::NeutralAtomMapper::findBestSwap() {
  auto swaps = getAllPossibleSwaps();
  // evaluate swaps based on cost function
  std::vector<std::pair<Swap, double>> swapCosts;
  for (const auto& swap : swaps) {
    swapCosts.emplace_back(swap, distanceCost(swap));
  }
  // get swap of minimal cost
  auto bestSwap = std::min_element(swapCosts.begin(), swapCosts.end(),
                                   [](const auto& swap1, const auto& swap2) {
                                     return swap1.second < swap2.second;
                                   });
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

fp NeutralAtomMapper::distanceCost(qc::Swap swap) {
  // compute the change in total distance
  fp distanceChange = 0;
  for (const auto& gate : this->frontLayer) {
    auto gateQubits = gate->get()->getUsedQubits();
    // distance before
    distanceChange -= this->hardwareQubits.getTotalDistance(gateQubits);

    // do swap
    if (gateQubits.find(swap.first) != gateQubits.end()) {
      gateQubits.insert(swap.second);
      gateQubits.erase(swap.first);
    } else if (gateQubits.find(swap.second) != gateQubits.end()) {
      gateQubits.insert(swap.first);
      gateQubits.erase(swap.second);
    } else {
      continue;
    }
    // distance after
    distanceChange += this->hardwareQubits.getTotalDistance(gateQubits);
  }
  return distanceChange;
}

} // namespace qc
