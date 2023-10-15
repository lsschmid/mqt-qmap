//
// Created by Ludwig Schmid on 05.10.23.
//

#pragma once

#include "NeutralAtomMappingResults.hpp"
#include "QuantumComputation.hpp"
#include "namap/NeutralAtomArchitecture.hpp"

namespace qc {
class NeutralAtomMapper {
protected:
  using Layer = std::set<std::unique_ptr<qc::Operation>*>;

  qc::NeutralAtomArchitecture arch;
  qc::QuantumComputation      mappedQc;
  Layer                       frontLayer;
  DAG                         dag;
  DAGIterators                dagIterators;
  //  qc::QuantumComputation      lookaheadLayer;
  //  std::map<Qubit, Qubit>      circToHwMapping;
  //  std::map<Qubit, CoordIndex> hwQubitCoordinates;
  NeutralAtomMappingResults results;

  struct Settings {
    fp a = 0.5;
  };

  void createFrontLayer();
  void
  updateFrontLayerByGate(std::set<std::unique_ptr<Operation>*>& gatesToExecute);
  void        updateFrontLayerByQubit(std::vector<Qubit>& qubitsToCheck);
  void        mapGate(std::unique_ptr<qc::Operation>* op);
  static bool isExecutable(qc::Operation* op);
  bool        isAtFront(std::unique_ptr<qc::Operation>* opPointer);

public:
  // Getter
  [[nodiscard]] inline qc::NeutralAtomArchitecture getArchitecture() const {
    return arch;
  }
  [[nodiscard]] inline NeutralAtomMappingResults getResults() const {
    return results;
  }
  [[nodiscard]] inline std::map<uint32_t, uint32_t> getInitialMapping() const {
    return mappedQc.initialLayout;
  }

  // Constructors
  NeutralAtomMapper(const qc::NeutralAtomArchitecture arch) : arch(arch){};

  // Methods
  void map(qc::QuantumComputation& qc);
};

} // namespace qc
