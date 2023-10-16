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
  std::vector<Qubit>          frontQubitsToUpdate;
  std::vector<Layer>          frontCandidates;
  DAG                         dag;
  DAGIterators                frontLayerIterators;
  Layer                       lookaheadLayer;
  std::vector<uint32_t>       lookaheadOffsets;
  std::vector<Layer>          lookaheadCandidates;
  std::vector<Qubit>          lookaheadQubitsToUpdate;
  uint32_t                    lookaheadDepth = 1;
  //  std::map<Qubit, Qubit>      circToHwMapping;
  //  std::map<Qubit, CoordIndex> hwQubitCoordinates;
  NeutralAtomMappingResults results;

  void createFrontLayer();
  void
  updateFrontLayerByGate(std::set<std::unique_ptr<Operation>*>& gatesToExecute);
  void        updateFrontLayerByQubit();
  void        updateFrontLayerByCandidates();
  void        findFrontCandidates();
  void        updateLookaheadLayerByQubit();
  void        findLookaheadCandidates();
  void        updateLookaheadLayerByCandidates();
  void        mapGate(std::unique_ptr<qc::Operation>* op);
  static bool commutesWith(Layer&                          layer,
                           std::unique_ptr<qc::Operation>* opPointer);
  static bool commute(std::unique_ptr<qc::Operation>* opPointer1,
                      std::unique_ptr<qc::Operation>* opPointer2);
  static bool
  commuteSingleQubitZAndControl(std::unique_ptr<qc::Operation>* opPointer1,
                                std::unique_ptr<qc::Operation>* opPointer2);
  static bool isExecutable(qc::Operation* op);
  void        addToFrontLayer(std::unique_ptr<qc::Operation>* opPointer);

  // temp
  void printLayers();

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
