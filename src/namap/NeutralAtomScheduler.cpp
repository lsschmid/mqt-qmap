//
// Created by Ludwig Schmid on 11.10.23.
//

#include "namap/NeutralAtomScheduler.hpp"

#include "CircuitOptimizer.hpp"

qc::SchedulerResults
qc::NeutralAtomScheduler::schedule(qc::QuantumComputation& qc, bool verbose) {
  CircuitOptimizer::decomposeSWAP(qc, false);
  if (verbose) {
    std::cout << "\n* schedule start!\n";
  }

  std::vector<fp> totalExecutionTimes(arch.getNpositions(), 0);
  std::vector<fp> blockedQubitsTimes(arch.getNpositions(), 0);
  fp              totalGateTime       = 0;
  fp              totalGateFidelities = 1;

  int index = 0;
  for (const auto& op : qc) {
    index++;
    if (verbose) {
      std::cout << "\n" << index << "\n";
    }

    auto qubits     = op->getUsedQubits();
    auto opTime     = arch.getOpTime(op.get());
    auto opFidelity = arch.getOpFidelity(op.get());

    // DEBUG info
    if (verbose) {
      std::cout << op->getName() << "  ";
      for (const auto& qubit : qubits) {
        std::cout << "q" << qubit << " ";
      }
      std::cout << "-> time: " << opTime << ", fidelity: " << opFidelity
                << "\n";
    }

    fp maxTime = 0;
    if (op->getType() != qc::AodMove && op->getType() != qc::AodActivate &&
        op->getType() != qc::AodDeactivate && qubits.size() > 1) {
      // multi qubit gates -> take into consideration blocking
      auto blockedQubits = arch.getBlockedCoordIndices(op.get());
      // get max execution time over all blocked qubits
      for (const auto& qubit : blockedQubits) {
        maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
        maxTime = std::max(maxTime, blockedQubitsTimes[qubit]);
      }

      // update total execution times
      for (const auto& qubit : qubits) {
        totalExecutionTimes[qubit] = maxTime + opTime;
      }
      for (const auto& qubit : blockedQubits) {
        blockedQubitsTimes[qubit] = maxTime + opTime;
      }

    } else {
      // other operations -> no blocking
      // get max execution time over all qubits
      for (const auto& qubit : qubits) {
        maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
      }
      // update total execution times
      for (const auto& qubit : qubits) {
        totalExecutionTimes[qubit] = maxTime + opTime;
      }
    }

    totalGateFidelities *= opFidelity;
    totalGateTime += opTime;
    if (verbose) {
      std::cout << "\n";
      printTotalExecutionTimes(totalExecutionTimes, blockedQubitsTimes);
    }
  }
  //  if (verbose) {
  std::cout << "\n* schedule end!\n";
  //}

  const auto maxExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  const auto totalIdleTime =
      maxExecutionTime * arch.getNpositions() - totalGateTime;
  const auto totalFidelities =
      totalGateFidelities *
      std::exp(-totalIdleTime / arch.getDecoherenceTime());

  printSchedulerResults(totalExecutionTimes, totalIdleTime, totalGateFidelities,
                        totalFidelities);

  return {maxExecutionTime, totalIdleTime, totalGateFidelities,
          totalFidelities};
}

void qc::NeutralAtomScheduler::printSchedulerResults(
    std::vector<fp>& totalExecutionTimes, fp totalIdleTime,
    fp totalGateFidelities, fp totalFidelities) {
  auto totalExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  std::cout << "\ntotalExecutionTimes: " << totalExecutionTime << "\n";
  std::cout << "totalIdleTime: " << totalIdleTime << "\n";
  std::cout << "totalGateFidelities: " << totalGateFidelities << "\n";
  std::cout << "totalFidelities: " << totalFidelities << "\n";
}

void qc::NeutralAtomScheduler::printTotalExecutionTimes(
    std::vector<fp>& totalExecutionTimes, std::vector<fp>& blockedQubitsTimes) {
  std::cout << "ExecutionTime: "
            << "\n";
  for (size_t qubit = 0; qubit < totalExecutionTimes.size(); qubit++) {
    std::cout << "[" << qubit << "] " << totalExecutionTimes[qubit] << " \t"
              << blockedQubitsTimes[qubit] << "\n";
  }
}
