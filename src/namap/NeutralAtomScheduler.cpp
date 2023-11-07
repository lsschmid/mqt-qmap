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

  int index = 0;
  for (auto& op : qc) {
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

    auto blockedQubits = arch.getBlockedCoordIndices(op.get());

    if (verbose) {
      std::cout << "blockedQubits: ";
      for (const auto& qubit : blockedQubits) {
        std::cout << "q" << qubit << " ";
      }
    }

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
    if (verbose) {
      std::cout << "\n";
      printTotalExecutionTimes(totalExecutionTimes);
    }
  }
  if (verbose) {
    printSchedulerResults(totalExecutionTimes, totalIdleTime, totalFidelities);
    std::cout << "\n* schedule end!\n";
  }

  return {totalExecutionTimes, totalIdleTime, totalFidelities};
}

void qc::NeutralAtomScheduler::printSchedulerResults(
    std::vector<fp>& totalExecutionTimes, fp totalIdleTime,
    fp totalFidelities) {
  auto totalExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  std::cout << "\ntotalExecutionTimes: " << totalExecutionTime << "\n";
  std::cout << "totalIdleTime: " << totalIdleTime << "\n";
  std::cout << "totalFidelities: " << totalFidelities << "\n";
}

void qc::NeutralAtomScheduler::printTotalExecutionTimes(
    std::vector<fp>& totalExecutionTimes) {
  std::cout << "ExecutionTime: "
            << "\n";
  int i = 0;
  for (auto qubit : totalExecutionTimes) {
    std::cout << "[" << i++ << "] " << qubit << "\n";
  }
}
