//
// Created by Ludwig Schmid on 05.10.23.
//

#include "namap/NeutralAtomMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "utils.hpp"

void qc::NeutralAtomMapper::map(qc::QuantumComputation& qc) {
  this->dag = qc::CircuitOptimizer::constructDAG(qc);
  for (auto& i : dag) {
    this->frontLayerIterators.push_back(i.begin());
    this->lookaheadOffsets.push_back(0);
    this->lookaheadCandidates.emplace_back();
    this->frontCandidates.emplace_back();
  }

  std::cout << "test" << std::endl;

  auto d = this->hardwareQubits.getSwapDistance(0, 1);
  this->hardwareQubits.move(0, 20, this->arch);
  d = this->hardwareQubits.getSwapDistance(0, 1);

  //   precompute all distances

  createFrontLayer();
  //  printLayers();
  updateLookaheadLayerByQubit();
  printLayers();

  // probably need
  // updateFrontLayerByQubit(list with all qubits) to trigger commutation

  // TODO: find SWAP
  // TODO: get list of front gates that can be executed now
  for (int i = 0; i < 6; i++) {
    std::set<std::unique_ptr<Operation>*> gatesToExecute;
    for (auto* gate : this->frontLayer) {
      gatesToExecute.insert(gate);
    }
    std::cout << i << std::endl;
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
      this->frontLayerIterators[gateQubit]++;
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
      auto* opPointer = *(this->frontLayerIterators[qubit]);
      auto* op        = opPointer->get();
      if (isExecutable(op)) {
        mapGate(opPointer);
        tempIter++;
        this->frontLayerIterators[qubit]++;
      } else {
        this->frontCandidates[qubit].insert(opPointer);
        tempIter++;

        // continue if following gates commute
        bool commutes = true;
        while (commutes && tempIter != this->dag[qubit].end()) {
          auto* nextOpPointer = *tempIter;
          auto* nextOp        = nextOpPointer->get();
          commutes            = commutesWith(this->frontLayer, nextOpPointer) &&
                     commutesWith(this->frontCandidates[qubit], nextOpPointer);
          if (commutes) {
            if (isExecutable(nextOp)) {
              mapGate(nextOpPointer);
              this->frontLayerIterators[qubit]++;
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
    Layer& layer, std::unique_ptr<qc::Operation>* opPointer) {
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
  std::cout << "mapped " << (*op)->getName();
  for (auto qubit : (*op)->getUsedQubits()) {
    std::cout << " " << qubit;
  }
  std::cout << std::endl;
}

bool qc::NeutralAtomMapper::isExecutable(qc::Operation* op) {
  if (op->getUsedQubits().size() == 1) {
    return true;
  }
  // TODO: check if two qubit operation is executable
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
