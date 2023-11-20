#include "namap/NeutralAtomScheduler.hpp"

#include "CircuitOptimizer.hpp"

qc::SchedulerResults
qc::NeutralAtomScheduler::schedule(qc::QuantumComputation& qc, bool verbose) {
  // decompose CX gates
  qc::CircuitOptimizer::replaceMCXWithMCZ(qc);

  if (verbose) {
    std::cout << "\n* schedule start!\n";
  }

  std::vector<fp> totalExecutionTimes(arch.getNpositions(), 0);
  // saves for each coord the time slots that are blocked by a multi qubit gate
  std::vector<std::deque<std::pair<fp, fp>>> blockedQubitsTimes(
      arch.getNpositions(), std::deque<std::pair<fp, fp>>());
  fp totalGateTime       = 0;
  fp totalGateFidelities = 1;

  int      index        = 0;
  int      nAodActivate = 0;
  uint32_t nCZs         = 0;
  for (const auto& op : qc) {
    index++;
    if (verbose) {
      std::cout << "\n" << index << "\n";
    }
    if (op->getType() == qc::AodActivate) {
      nAodActivate++;
    } else if (op->getType() == qc::OpType::Z && op->getNcontrols() == 1) {
      nCZs++;
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
      bool blocked = true;
      while (blocked) {
        // get regular max execution time
        for (const auto& qubit : qubits) {
          maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
        }
        // check if all blocked qubits are free at maxTime
        blocked = false;
        for (const auto& qubit : blockedQubits) {
          // check if qubit is blocked at maxTime
          if (!blockedQubitsTimes[qubit].empty() &&
              blockedQubitsTimes[qubit].front().first < maxTime &&
              blockedQubitsTimes[qubit].front().second > maxTime) {
            blocked = true;
            // update maxTime to the end of the blocking
            maxTime = blockedQubitsTimes[qubit].front().second;
            // remove the blocking
            break;
          }
        }
      }

      // update total execution times
      for (const auto& qubit : qubits) {
        totalExecutionTimes[qubit] = maxTime + opTime;
      }
      for (const auto& qubit : blockedQubits) {
        blockedQubitsTimes[qubit].emplace_back(maxTime, maxTime + opTime);
      }

    } else {
      // other operations -> no blocking
      // get max execution time over all qubits
      for (const auto& qubit : qubits) {
        maxTime = std::max(maxTime, totalExecutionTimes[qubit]);
        // remove all blocked times that are smaller than maxTime
        while (!blockedQubitsTimes[qubit].empty() &&
               blockedQubitsTimes[qubit].front().second < maxTime) {
          blockedQubitsTimes[qubit].pop_front();
        }
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
  if (verbose) {
    std::cout << "nAodActivate: " << nAodActivate << "\n";
  }

  const auto maxExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  const auto totalIdleTime =
      maxExecutionTime * arch.getNqubits() - totalGateTime;
  const auto totalFidelities =
      totalGateFidelities *
      std::exp(-totalIdleTime / arch.getDecoherenceTime());

  printSchedulerResults(totalExecutionTimes, totalIdleTime, totalGateFidelities,
                        totalFidelities, nCZs);

  return {maxExecutionTime, totalIdleTime, totalGateFidelities, totalFidelities,
          nCZs};
}

void qc::NeutralAtomScheduler::printSchedulerResults(
    std::vector<fp>& totalExecutionTimes, fp totalIdleTime,
    fp totalGateFidelities, fp totalFidelities, uint32_t nCZs) {
  auto totalExecutionTime =
      *std::max_element(totalExecutionTimes.begin(), totalExecutionTimes.end());
  std::cout << "\ntotalExecutionTimes: " << totalExecutionTime << "\n";
  std::cout << "totalIdleTime: " << totalIdleTime << "\n";
  std::cout << "totalGateFidelities: " << totalGateFidelities << "\n";
  std::cout << "totalFidelities: " << totalFidelities << "\n";
  std::cout << "totalnCZs: " << nCZs << "\n";
}

void qc::NeutralAtomScheduler::printTotalExecutionTimes(
    std::vector<fp>&                            totalExecutionTimes,
    std::vector<std::deque<std::pair<fp, fp>>>& blockedQubitsTimes) {
  std::cout << "ExecutionTime: "
            << "\n";
  for (size_t qubit = 0; qubit < totalExecutionTimes.size(); qubit++) {
    std::cout << "[" << qubit << "] " << totalExecutionTimes[qubit] << " \t";
    for (const auto& blockedTime : blockedQubitsTimes[qubit]) {
      std::cout << blockedTime.first << "-" << blockedTime.second << " \t";
    }
    std::cout << "\n";
  }
}
