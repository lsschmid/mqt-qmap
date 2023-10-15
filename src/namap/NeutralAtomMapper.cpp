//
// Created by Ludwig Schmid on 05.10.23.
//

#include "namap/NeutralAtomMapper.hpp"

#include "CircuitOptimizer.hpp"
#include "utils.hpp"

void qc::NeutralAtomMapper::map(qc::QuantumComputation& qc) {
  this->dag = qc::CircuitOptimizer::constructDAG(qc);
  for (auto& i : dag) {
    this->dagIterators.push_back(i.begin());
  }
  this->mappedQc = qc::QuantumComputation(qc.getNqubits());

  createFrontLayer();

  // TODO: find SWAP
  // TODO: get list of front gates that can be executed now
  for (int i = 0; i < 6; i++) {
    std::set<std::unique_ptr<Operation>*> gatesToExecute;
    for (auto* gate : this->frontLayer) {
      gatesToExecute.insert(gate);
    }

    updateFrontLayerByGate(gatesToExecute);
  }
}

void qc::NeutralAtomMapper::createFrontLayer() {
  // vector indicating if qubit is free
  std::vector<Qubit> qubitsToCheck;
  for (uint32_t i = 0; i < this->dag.size(); i++) {
    qubitsToCheck.push_back(i);
  }
  updateFrontLayerByQubit(qubitsToCheck);
}

void qc::NeutralAtomMapper::updateFrontLayerByGate(
    std::set<std::unique_ptr<Operation>*>& gatesToExecute) {
  std::vector<Qubit> qubitsToCheck;
  for (auto* gate : gatesToExecute) {
    // add other gate qubits to be checked
    for (auto qubit : (*gate)->getUsedQubits()) {
      if (std::find(qubitsToCheck.begin(), qubitsToCheck.end(), qubit) ==
          qubitsToCheck.end()) {
        qubitsToCheck.push_back(qubit);
      }
    }
    // remove from current FrontLayer
    mapGate(gate);
    // update DAG iterators
    for (auto gateQubit : (*gate)->getUsedQubits()) {
      this->dagIterators[gateQubit]++;
    }
    // remove from FrontLayer
    this->frontLayer.erase(
        std::find(this->frontLayer.begin(), this->frontLayer.end(), gate));
  }
  updateFrontLayerByQubit(qubitsToCheck);
}

void qc::NeutralAtomMapper::updateFrontLayerByQubit(
    std::vector<Qubit>& qubitsToCheck) {
  for (size_t i = 0; i < qubitsToCheck.size(); i++) {
    auto qubit = qubitsToCheck[i];
    while (this->dagIterators[qubit] != this->dag[qubit].end()) {
      auto* opPointer = *(this->dagIterators[qubit]);
      auto* op        = opPointer->get();
      // add to mapped qc if executable
      if (isExecutable(op)) {
        mapGate(opPointer);
        this->dagIterators[qubit]++;
      } else if (isAtFront(opPointer)) {
        // add to front layer if not yet there
        this->frontLayer.insert(opPointer);

        break;
      } else {
        // add other gate qubits to be checked
        for (auto gateQubit : op->getUsedQubits()) {
          if (gateQubit != qubit &&
              std::find(qubitsToCheck.begin(), qubitsToCheck.end(),
                        gateQubit) == qubitsToCheck.end()) {
            qubitsToCheck.push_back(gateQubit);
          }
        }
        break;
      }
    }
  }
}

void qc::NeutralAtomMapper::mapGate(std::unique_ptr<qc::Operation>* op) {
  this->mappedQc.emplace_back((*op)->clone());
}

bool qc::NeutralAtomMapper::isExecutable(qc::Operation* op) {
  if (op->getUsedQubits().size() == 1) {
    return true;
  }
  // TODO: check if two qubit operation is executable
  return false;
}

bool qc::NeutralAtomMapper::isAtFront(
    std::unique_ptr<qc::Operation>* opPointer) {
  auto* op = opPointer->get();
  for (auto qubit : op->getUsedQubits()) {
    if (this->dagIterators[qubit] != std::find(this->dag[qubit].begin(),
                                               this->dag[qubit].end(),
                                               opPointer)) {
      return false;
    }
  }
  return true;
}
