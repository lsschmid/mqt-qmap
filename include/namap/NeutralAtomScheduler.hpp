//
// Created by Ludwig Schmid on 11.10.23.
//
#pragma once

#include "QuantumComputation.hpp"
#include "namap/HardwareQubits.hpp"
#include "namap/NeutralAtomArchitecture.hpp"

namespace qc {
/**
 * @brief Struct to store the results of the scheduler
 */
struct SchedulerResults {
  fp       totalExecutionTime;
  fp       totalIdleTime;
  fp       totalGateFidelities;
  fp       totalFidelities;
  uint32_t nCZs = 0;

  SchedulerResults(fp totalExecutionTime, fp totalIdleTime,
                   fp totalGateFidelities, fp totalFidelities, uint32_t nCZs)
      : totalExecutionTime(totalExecutionTime), totalIdleTime(totalIdleTime),
        totalGateFidelities(totalGateFidelities),
        totalFidelities(totalFidelities), nCZs(nCZs) {}

  std::string inline toString() {
    std::stringstream ss;
    ss << "Total execution time: " << totalExecutionTime;
    ss << "\nTotal idle time: " << totalIdleTime
       << "\nTotal fidelities: " << totalFidelities;
    return ss.str();
  }
  std::string inline toCsv() {
    std::stringstream ss;
    ss << totalExecutionTime << ", " << totalIdleTime << "," << totalFidelities;
    return ss.str();
  }
};

/**
 * @brief Class to schedule a quantum circuit on a neutral atom architecture
 * @details For each gate/operation in the input circuit, the scheduler checks the
 * earliest possible time slot for execution. If the gate is a multi qubit gate,
 * also the blocking of other qubits is taken into consideration.
 * The execution times are read from the neutral atom architecture.
 */
class NeutralAtomScheduler {
protected:
  qc::NeutralAtomArchitecture arch;

public:
  // Constructor
  NeutralAtomScheduler(const qc::NeutralAtomArchitecture& arch) : arch(arch) {}

  /**
   * @brief Schedules the given quantum circuit on the neutral atom architecture
   * @details For each gate/operation in the input circuit, the scheduler checks the
   * earliest possible time slot for execution. If the gate is a multi qubit gate,
   * also the blocking of other qubits is taken into consideration.
   * The execution times are read from the neutral atom architecture.
   * @param qc Quantum circuit to schedule
   * @param verbose If true, prints additional information
   * @return SchedulerResults
   */
  SchedulerResults schedule(qc::QuantumComputation& qc, bool verbose);

  // Helper Print functions
  static void printSchedulerResults(std::vector<fp>& totalExecutionTimes,
                                    fp totalIdleTime, fp totalGateFidelities,
                                    fp totalFidelities, uint32_t nCZs);
  static void printTotalExecutionTimes(
      std::vector<fp>&                            totalExectuionTimes,
      std::vector<std::deque<std::pair<fp, fp>>>& blockedQubitsTimes);
};

} // namespace qc
