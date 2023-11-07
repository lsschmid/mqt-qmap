//
// Created by Ludwig Schmid on 11.10.23.
//
#pragma once

#include "QuantumComputation.hpp"
#include "namap/HardwareQubits.hpp"
#include "namap/NeutralAtomArchitecture.hpp"

namespace qc {
struct SchedulerResults {
  std::vector<fp> totalExecutionTimes;
  fp              totalIdleTime;
  fp              totalFidelities;

  SchedulerResults(std::vector<fp> totalExecutionTimes, fp totalIdleTime,
                   fp totalFidelities)
      : totalExecutionTimes(std::move(totalExecutionTimes)),
        totalIdleTime(totalIdleTime), totalFidelities(totalFidelities) {}

  std::string inline toString() {
    std::stringstream ss;
    ss << "Total execution times: ";
    for (auto& time : totalExecutionTimes) {
      ss << time << " ";
    }
    ss << "\nTotal idle time: " << totalIdleTime
       << "\nTotal fidelities: " << totalFidelities;
    return ss.str();
  }
  std::string inline toCsv() {
    std::stringstream ss;
    for (auto& time : totalExecutionTimes) {
      ss << time << ",";
    }
    ss << totalIdleTime << "," << totalFidelities;
    return ss.str();
  }
};

class NeutralAtomScheduler {
protected:
  qc::NeutralAtomArchitecture arch;
  std::vector<fp>             totalExecutionTimes;
  fp                          totalIdleTime;
  fp                          totalFidelities;

public:
  NeutralAtomScheduler(const qc::NeutralAtomArchitecture& arch)
      : arch(arch), totalExecutionTimes(std::vector<fp>(arch.getNqubits(), 0)),
        totalIdleTime(0), totalFidelities(1.0){};

  SchedulerResults schedule(qc::QuantumComputation& qc, bool verbose);

  static void printSchedulerResults(std::vector<fp>& totalExecutionTimes,
                                    fp totalIdleTime, fp totalFidelities);
  static void printTotalExecutionTimes(std::vector<fp>& totalExectuionTimes);
  static void printQubitsInfo(std::set<unsigned int> qubits,
                              qc::HardwareQubits     hardwareQubits);
};

} // namespace qc
