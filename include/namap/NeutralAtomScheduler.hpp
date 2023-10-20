//
// Created by Ludwig Schmid on 11.10.23.
//
#pragma once

#include "QuantumComputation.hpp"
#include "namap/NeutralAtomArchitecture.hpp"
#include "namap/NeutralAtomMapper.hpp"
#include "namap/NeutralAtomUtils.hpp"

namespace qc {

class NeutralAtomScheduler {
protected:
  qc::NeutralAtomArchitecture arch;
  qc::HardwareQubits          hardwareQubits;
  InitialMapping              mappingType;
  qc::Mapping                 mapping;
  std::vector<fp>             totalExecutionTimes;
  fp                          totalIdleTime;
  fp                          totalFidelities;

  // Methods

  // temp`

public:
  // Getter
  [[nodiscard]] inline qc::NeutralAtomArchitecture getArchitecture() const {
    return arch;
  }

  // Constructors
  NeutralAtomScheduler(const qc::NeutralAtomArchitecture& arch,
                       qc::QuantumComputation& qc, InitialMapping mappingType,
                       InitialCoordinateMapping initialCoordinateMapping)
      : arch(arch), hardwareQubits(arch, initialCoordinateMapping),
        mappingType(mappingType), mapping(arch.getNqubits(), mappingType),
        totalExecutionTimes(std::vector<fp>(qc.getNqubits(), 0)),
        totalIdleTime(0), totalFidelities(1.0){};

  // Methods
  void schedule(qc::QuantumComputation& qc, bool isBlockedForAll);
};

} // namespace qc
