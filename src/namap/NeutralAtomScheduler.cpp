//
// Created by Ludwig Schmid on 11.10.23.
//

#include "namap/NeutralAtomScheduler.hpp"

void qc::NeutralAtomScheduler::schedule(qc::QuantumComputation& qc,
                                        bool isBlockedForAll) {
  std::cout << "\n* schedule start!" << std::endl;
  // search the gates of the circuit in order
  int index = 0;
  for (auto& op : qc) {
    double opTime, opFidelity;

    // calculate time & fidelity
    // TODO: make "shuttling" operation in the .QASM
    auto qubits = op->getUsedQubits();

    std::cout << "\n" << index++ << std::endl;
    if (op->getType() != None) {
      // gate
      std::cout << op->getName() << "  ";
      opTime     = arch.getGateTime(op->getType());
      opFidelity = arch.getGateAverageFidelity(op->getType());

      // circuitQubit -> hardwareQubit -> coordinates
      for (auto& qubit : qubits) {
        // circuit qubits
        std::cout << "q" << qubit << "{";
        // hardware qubits
        std::cout << mapping.getHwQubit(qubit) << "(";
        // coordinate index
        std::cout << hardwareQubits.getCoordIndex(mapping.getHwQubit(qubit))
                  << ")}, ";
      }
      std::cout << std::endl;
      std::cout << "-> time: " << opTime << ", fidelity: " << opFidelity
                << std::endl;
    } else {
      // shuttling
      /*
      opTime     = arch.getShuttlingTime(
          static_cast<ShuttlingType>(op->getShuttlingType()));
      opFidelity = arch.getShuttlingAverageFidelity(
          static_cast<ShuttlingType>(op->getShuttlingType()));
      */
    }

    // TODO: calculate the number of circuit qubits
    // op->getBlockedQubits;
    std::set<unsigned int> blockedQubits = qubits;
    if (op->getType() != None && qubits.size() != 1) {
      for (auto q = 0; q < qc.getNqubits(); q++) {
        auto qHardwareQubit = mapping.getHwQubit(q);
        auto qCoordIndex    = hardwareQubits.getCoordIndex(qHardwareQubit);
        for (auto& qubit : qubits) {
          if (q == qubit) {
            continue;
          }
          auto qubitHardwareQubit = mapping.getHwQubit(qubit);
          auto qubitCoordIndex =
              hardwareQubits.getCoordIndex(qubitHardwareQubit);
          auto qubitDistance =
              arch.getEuclidianDistance(qCoordIndex, qubitCoordIndex);
          if (qubitDistance < arch.getInteractionRadius()) {
            blockedQubits.insert(q);
          }
        }
      }
    } else if (op->getType() == None) {
      // TODO: getBlockedQubits for Shuttling operations
      if (!isBlockedForAll) {
        // calculate row and col
        for (auto& qubit : qubits) {
          auto qHardwareQubit = mapping.getHwQubit(qubit);
          auto qCoordIndex    = hardwareQubits.getCoordIndex(qHardwareQubit);
          auto qCoord         = arch.getCoordinate(qCoordIndex).getXY();

          // for all Rows
          for (int i = 0; i < arch.getNcolumns(); i++) {
            blockedQubits.insert(qCoord.first * arch.getNcolumns() + i);
          }
          // for all Columns
          for (int i = 0; i < arch.getNrows(); i++) {
            blockedQubits.insert(i * arch.getNcolumns() + qCoord.second);
          }
        }
        // Move vs. Activate/Deactivate

      } else {
        for (auto i = 0; i < qc.getNqubits(); i++) {
          blockedQubits.insert(i);
        }
      }
    }
    std::cout << "blockedQubits: ";
    for (auto& qubit : blockedQubits) {
      std::cout << qubit << ", ";
    }
    std::cout << std::endl;

    // find maxTime for blockedQubits
    auto maxQubit = *std::max_element(
        blockedQubits.begin(), blockedQubits.end(), [this](int a, int b) {
          return totalExecutionTimes[a] < totalExecutionTimes[b];
        });
    fp maxTime = totalExecutionTimes[maxQubit];

    // update totalExecutionTimes & totalIdleTime
    for (auto qubit : blockedQubits) {
      if (totalExecutionTimes[qubit] < maxTime) {
        totalIdleTime += maxTime - totalExecutionTimes[qubit];
      }
      totalExecutionTimes[qubit] = maxTime + opTime;
    }
    std::cout << "ExecutionTime: " << std::endl;
    int i = 0;
    for (auto qubit : totalExecutionTimes) {
      std::cout << "[" << i << "] " << qubit << std::endl;
      i++;
    }

    // update totalFidelities
    totalFidelities *= opFidelity;
  }

  // print the scheduler's results
  auto maxIter =
      std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  fp totalExecutionTime = *maxIter;
  std::cout << "\ntotalExecutionTimes: " << totalExecutionTime << std::endl;
  std::cout << "totalIdleTime: " << totalIdleTime << std::endl;
  std::cout << "totalFidelities: " << totalFidelities << std::endl;
}
