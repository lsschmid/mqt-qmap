//
// Created by Ludwig Schmid on 11.10.23.
//
#pragma once

#include "QuantumComputation.hpp"
#include "namap/HardwareQubits.hpp"
#include "namap/NeutralAtomArchitecture.hpp"

namespace qc {

class NeutralAtomScheduler {
protected:
  qc::NeutralAtomArchitecture arch;
  qc::HardwareQubits          hardwareQubits;
  std::vector<fp>             totalExecutionTimes;
  fp                          totalIdleTime;
  fp                          totalFidelities;

public:
  [[nodiscard]] inline qc::NeutralAtomArchitecture getArchitecture() const {
    return arch;
  }

  NeutralAtomScheduler(const qc::NeutralAtomArchitecture& arch,
                       InitialCoordinateMapping initialCoordinateMapping)
      : arch(arch), hardwareQubits(arch, initialCoordinateMapping),
        totalExecutionTimes(std::vector<fp>(arch.getNqubits(), 0)),
        totalIdleTime(0), totalFidelities(1.0){};

  void schedule(qc::QuantumComputation& qc, bool isBlockedForAll);

  static void printSchedulerResults(std::vector<fp>& totalExecutionTimes,
                                    fp totalIdleTime, fp totalFidelities);
  static void printTotalExecutionTimes(std::vector<fp>& totalExectuionTimes);
  static void printQubitsInfo(std::set<unsigned int> qubits,
                              qc::HardwareQubits     hardwareQubits);
};

} // namespace qc
