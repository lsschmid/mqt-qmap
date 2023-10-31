//
// Created by Ludwig Schmid on 05.10.23.
//

#pragma once

#include "NeutralAtomMappingResults.hpp"
#include "QuantumComputation.hpp"
#include "namap/AodScheduler.hpp"
#include "namap/HardwareQubits.hpp"
#include "namap/Mapping.hpp"
#include "namap/NeutralAtomArchitecture.hpp"
#include "namap/NeutralAtomUtils.hpp"

namespace qc {

using GateList = std::set<std::unique_ptr<qc::Operation>*>;
class NeutralAtomMapper {
protected:
  struct MapperParameters {
    fp lookaheadWeightSwaps         = 0.1;
    fp lookaheadWeightMoves         = 0.1;
    fp decay                        = 0.1;
    fp shuttlingTimeWeight          = 1;
    fp shuttlingMakeExecutableBonus = 1;
  };

  qc::NeutralAtomArchitecture                  arch;
  qc::QuantumComputation                       mappedQc;
  std::vector<std::unique_ptr<qc::Operation>*> executedCommutingGates;
  GateList                                     frontLayer;
  std::set<Qubit>                              frontQubitsToUpdate;
  std::vector<GateList>                        frontCandidates;
  DAG                                          dag;
  DAGIterators                                 frontLayerIterators;
  GateList                                     lookaheadLayer;
  std::vector<uint32_t>                        lookaheadOffsets;
  std::vector<GateList>                        lookaheadCandidates;
  std::set<Qubit>                              lookaheadQubitsToUpdate;
  uint32_t                                     lookaheadDepth = 1;
  MapperParameters                             parameters;
  std::deque<std::set<HwQubit>>                lastBlockedQubits;
  std::deque<AtomMove>                         lastMoves;
  std::vector<fp>                              decayWeights;
  uint32_t                                     nSwaps = 0;
  uint32_t                                     nMoves = 0;

  //  NeutralAtomMappingResults   results;
  HardwareQubits hardwareQubits;
  Mapping        mapping;
  bool           verbose = true;

  // Methods for layer creation
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
  static bool commutesWith(GateList&                       layer,
                           std::unique_ptr<qc::Operation>* opPointer);
  static bool commute(std::unique_ptr<qc::Operation>* opPointer1,
                      std::unique_ptr<qc::Operation>* opPointer2);
  static bool
       commuteSingleQubitZAndControl(std::unique_ptr<qc::Operation>* opPointer1,
                                     std::unique_ptr<qc::Operation>* opPointer2);
  bool isExecutable(std::unique_ptr<qc::Operation>* opPointer);
  void addToFrontLayer(std::unique_ptr<qc::Operation>* opPointer);

  // Methods for mapping
  Swap               findBestSwap();
  GateList           getExecutableGates();
  std::set<Swap>     getAllPossibleSwaps();
  AtomMove           findBestAtomMove();
  std::set<MoveComb> getAllPossibleMoveCombinations();
  std::set<MoveComb> getNearbyMoveCombinations(HwQubit start, HwQubit target);
  std::set<MoveComb> getMoveAwayCombinations(HwQubit start, HwQubit target);

  void updateMapping(Swap swap);
  void updateMappingMove(AtomMove move);

  // Methods cost functions
  fp distanceCost(const Swap& swap);
  fp distancePerLayer(const Swap& swap, GateList& layer);
  fp moveCost(const AtomMove& move);
  fp moveCostComb(const MoveComb& moveComb);
  fp moveDistancePerLayer(const AtomMove& move, GateList& layer);
  fp parallelMoveCost(const AtomMove& move);

  // temp
  void printLayers();

public:
  // Getter
  [[nodiscard]] inline qc::NeutralAtomArchitecture getArchitecture() const {
    return arch;
  }
  //  [[nodiscard]] inline NeutralAtomMappingResults getResults() const {
  //    return results;
  //  }
  [[nodiscard]] inline std::map<uint32_t, uint32_t> getInitialMapping() const {
    return mappedQc.initialLayout;
  }

  // Constructors
  NeutralAtomMapper(const qc::NeutralAtomArchitecture& arch,
                    InitialCoordinateMapping           initialCoordinateMapping)
      : arch(arch), mappedQc(arch.getNpositions()),
        hardwareQubits(arch, initialCoordinateMapping){};

  // Methods
  QuantumComputation map(qc::QuantumComputation& qc);
};

} // namespace qc
