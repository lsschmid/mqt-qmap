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
