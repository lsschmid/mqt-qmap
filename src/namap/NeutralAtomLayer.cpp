//
// Created by Ludwig Schmid on 21.01.24.
//

#pragma once

#include "namap/NeutralAtomLayer.hpp"

#include "Definitions.hpp"
#include "namap/NeutralAtomDefinitions.hpp"

#include <algorithm>
#include <set>

namespace qc {

void NeutralAtomLayer::updateByQubits(const std::set<Qubit>& qubitsToUpdate) {
  updateCandidatesByQubits(qubitsToUpdate);
  candidatesToGates(qubitsToUpdate);
}

std::vector<uint32_t> NeutralAtomLayer::getIteratorOffset() {
  std::vector<uint32_t> offset;
  for (uint32_t i = 0; i < this->dag.size(); ++i) {
    offset.push_back(static_cast<uint32_t>(
        std::distance(this->dag[i].begin(), this->iterators[i])));
  }
  return offset;
}

void NeutralAtomLayer::initLayerOffset(
    const std::vector<uint32_t>& iteratorOffset) {
  this->gates.clear();
  for (auto& candidate : this->candidates) {
    candidate.clear();
  }
  this->MappedSingleQubitGates.clear();
  // if iteratorOffset is empty, set all iterators to begin
  if (iteratorOffset.empty()) {
    for (uint32_t i = 0; i < this->dag.size(); ++i) {
      this->iterators[i] = this->dag[i].begin();
    }
  } else {
    for (uint32_t i = 0; i < this->dag.size(); ++i) {
      this->iterators[i] = this->dag[i].begin() + iteratorOffset[i];
    }
  }
  std::set<Qubit> allQubits;
  for (uint32_t i = 0; i < this->dag.size(); ++i) {
    allQubits.insert(i);
  }
  updateByQubits(allQubits);
}

void NeutralAtomLayer::updateCandidatesByQubits(
    const std::set<Qubit>& qubitsToUpdate) {
  for (const auto& qubit : qubitsToUpdate) {
    auto tempIter = iterators[qubit];
    while (tempIter < this->dag[qubit].end()) {
      auto* op = (*tempIter)->get();
      if (op->getUsedQubits().size() == 1) {
        MappedSingleQubitGates.push_back(op);
        this->iterators[qubit]++;
        tempIter++;
      } else {
        // continue if following gates commute
        bool commutes = true;
        while (commutes && tempIter < this->dag[qubit].end()) {
          auto* nextOp = (*tempIter)->get();
          commutes     = commutesWithAtQubit(gates, nextOp, qubit) &&
                     commutesWithAtQubit(candidates[qubit], nextOp, qubit);
          if (commutes) {
            if (nextOp->getUsedQubits().size() == 1) {
              MappedSingleQubitGates.push_back(nextOp);
            } else { // not executable but commutes
              candidates[qubit].push_back(nextOp);
            }
          }
          tempIter++;
        }
        break;
      }
    }
  }
}

void NeutralAtomLayer::candidatesToGates(
    const std::set<Qubit>& qubitsToUpdate) {
  for (const auto& qubit : qubitsToUpdate) {
    std::vector<const Operation*> toRemove;
    for (const auto* opPointer : candidates[qubit]) {
      // check if gate is candidate for all qubits it uses
      bool inFrontLayer = true;
      for (const auto& opQubit : opPointer->getUsedQubits()) {
        if (qubit == opQubit) {
          continue;
        }
        if (std::find(candidates[opQubit].begin(), candidates[opQubit].end(),
                      opPointer) == candidates[opQubit].end()) {
          inFrontLayer = false;
          break;
        }
      }
      if (inFrontLayer) {
        this->gates.push_back(opPointer);
        // remove from candidacy of other qubits
        for (const auto& opQubit : opPointer->getUsedQubits()) {
          if (qubit == opQubit) {
            continue;
          }
          candidates[opQubit].erase(std::find(candidates[opQubit].begin(),
                                              candidates[opQubit].end(),
                                              opPointer));
        }

        // save to remove from candidacy of this qubit
        toRemove.push_back(opPointer);
      }
    }
    // remove from candidacy of this qubit
    // has to be done now to not change iterating list
    for (const auto* opPointer : toRemove) {
      candidates[qubit].erase(std::find(candidates[qubit].begin(),
                                        candidates[qubit].end(), opPointer));
    }
  }
}

void qc::NeutralAtomLayer::removeGatesAndUpdate(const GateList& gatesToRemove) {
  this->MappedSingleQubitGates.clear();
  std::set<Qubit> qubitsToUpdate;
  for (const auto& gate : gatesToRemove) {
    if (std::find(gates.begin(), gates.end(), gate) != gates.end()) {
      gates.erase(std::find(gates.begin(), gates.end(), gate));
      auto usedQubits = gate->getUsedQubits();
      qubitsToUpdate.insert(usedQubits.begin(), usedQubits.end());
    }
  }
  for (const auto& qubit : qubitsToUpdate) {
    ++this->iterators[qubit];
  }
  updateByQubits(qubitsToUpdate);
}

// Commutation

bool qc::NeutralAtomLayer::commutesWithAtQubit(const qc::GateList& layer,
                                               const Operation*    opPointer,
                                               const Qubit&        qubit) {
  return std::all_of(layer.begin(), layer.end(),
                     [&opPointer, &qubit](const auto& frontOpPointer) {
                       return commuteAtQubit(opPointer, frontOpPointer, qubit);
                     });
}

bool qc::NeutralAtomLayer::commuteAtQubit(const Operation* op1,
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

} // namespace qc