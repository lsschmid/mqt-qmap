//
// Created by Ludwig Schmid on 11.10.23.
//

#include "namap/NeutralAtomScheduler.hpp"

void qc::NeutralAtomScheduler::schedule(qc::QuantumComputation& qc,
                                        bool isBlockedForAll) {
  std::cout << "\n* schedule start!" << std::endl;

  int index = 0;
  for (auto& op : qc) {
    std::cout << "\n" << index++ << std::endl;

    auto   qubits = op->getUsedQubits();
    double opTime, opFidelity;

    // calculate time & fidelity
    auto opType = op->getType();
    if (opType == OpType::Move || opType == OpType::AodActivate ||
        opType == OpType::AodDeactivate) { // shuttling
      opTime     = arch.getShuttlingTime(op->getType());
      opFidelity = arch.getShuttlingAverageFidelity(op->getType());
    } else { // gate
      opTime     = arch.getGateTime(op->getType());
      opFidelity = arch.getGateAverageFidelity(op->getType());
    }

    // DEBUG info
    std::cout << op->getName() << "  ";
    printQubitsInfo(qubits, hardwareQubits);
    std::cout << "-> time: " << opTime << ", fidelity: " << opFidelity
              << std::endl;

    // getBlockedQubits;
    std::set<unsigned int> blockedQubits;
    if (qubits.size() == 1) {
      blockedQubits.insert(*qubits.begin());
    } else {
      blockedQubits = hardwareQubits.getBlockedQubits(qubits);
    }
    if (isBlockedForAll) {
      for (uint32_t i = 0; i < arch.getNqubits(); i++) {
        blockedQubits.insert(i);
      }
    }

    /*
    std::set<unsigned int> blockedQubits = qubits;
    if(!arch.checkShuttlingType(op->getType()) && qubits.size()!=1){
      for(auto q=0; q<arch.getNqubits(); q++){
        auto qCoordIndex    = hardwareQubits.getCoordIndex(q);
        for(auto& qubit : qubits){
          if(q==qubit) continue;

          auto qubitCoordIndex    = hardwareQubits.getCoordIndex(qubit);
          auto qubitDistance      = arch.getEuclidianDistance(qCoordIndex,
    qubitCoordIndex); if(qubitDistance < arch.getInteractionRadius()){
            blockedQubits.insert(q);
          }
        }
      }
    }
    else if(arch.checkShuttlingType(op->getType()) && isBlockedForAll){
      for(auto i=0; i<arch.getNqubits(); i++){
        blockedQubits.insert(i);
      }
    }
     */
    std::cout << "blockedQubits: ";
    printQubitsInfo(blockedQubits, hardwareQubits);

    // update totalExecutionTimes & totalIdleTime & totalFidelities
    auto maxQubit = *std::max_element(
        blockedQubits.begin(), blockedQubits.end(), [this](size_t a, size_t b) {
          return totalExecutionTimes[a] < totalExecutionTimes[b];
        });
    fp const maxTime = totalExecutionTimes[maxQubit];

    for (auto qubit : blockedQubits) {
      totalIdleTime += std::max(0.0, maxTime - totalExecutionTimes[qubit]);
      totalExecutionTimes[qubit] = maxTime + opTime;
    }
    totalFidelities *= opFidelity;
    printTotalExecutionTimes(totalExecutionTimes);
  }
  printSchedulerResults(totalExecutionTimes, totalIdleTime, totalFidelities);
}

void qc::NeutralAtomScheduler::printSchedulerResults(
    std::vector<fp>& totalExecutionTimes, fp totalIdleTime,
    fp totalFidelities) {
  auto totalExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  std::cout << "\ntotalExecutionTimes: " << totalExecutionTime << std::endl;
  std::cout << "totalIdleTime: " << totalIdleTime << std::endl;
  std::cout << "totalFidelities: " << totalFidelities << std::endl;
}

void qc::NeutralAtomScheduler::printTotalExecutionTimes(
    std::vector<fp>& totalExecutionTimes) {
  std::cout << "ExecutionTime: " << std::endl;
  int i = 0;
  for (auto qubit : totalExecutionTimes) {
    std::cout << "[" << i++ << "] " << qubit << std::endl;
  }
}

void qc::NeutralAtomScheduler::printQubitsInfo(
    std::set<unsigned int> qubits, qc::HardwareQubits hardwareQubits) {
  for (auto& qubit : qubits) {
    std::cout << "q" << qubit << "{";
    std::cout << hardwareQubits.getCoordIndex(qubit) << "}, ";
  }
  std::cout << std::endl;
}
