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

using GateList = std::vector<const Operation*>;

struct MapperParameters {
  fp lookaheadWeightSwaps = 0.1;
  fp lookaheadWeightMoves = 0.1;
  fp decay                = 0.1;
  fp shuttlingTimeWeight  = 1;
  fp gateWeight           = 1;
  fp shuttlingWeight      = 1;
};

struct MultiQubitMovePos {
  CoordIndices coords;
  size_t       nMoves;

  // copy constructor
  MultiQubitMovePos() : coords(), nMoves(0) {}
  MultiQubitMovePos(const MultiQubitMovePos& other)
      : coords(other.coords), nMoves(other.nMoves) {}

  // copy assignment
  MultiQubitMovePos& operator=(const MultiQubitMovePos& other) {
    coords = other.coords;
    nMoves = other.nMoves;
    return *this;
  }
};

class NeutralAtomMapper {
protected:
  qc::NeutralAtomArchitecture            arch;
  qc::QuantumComputation                 mappedQc;
  std::vector<const Operation*>          executedCommutingGates;
  GateList                               frontLayerGate;
  GateList                               frontLayerShuttling;
  std::set<Qubit>                        frontQubitsToUpdate;
  std::vector<GateList>                  frontCandidates;
  DAG                                    dag;
  DAGIterators                           frontLayerIterators;
  std::vector<SwapOrMove>                swapCloseByFront;
  std::vector<std::pair<SwapOrMove, fp>> moveExactFront;
  GateList                               lookaheadLayerGate;
  GateList                               lookaheadLayerShuttling;
  std::vector<uint32_t>                  lookaheadOffsets;
  std::vector<GateList>                  lookaheadCandidates;
  std::set<Qubit>                        lookaheadQubitsToUpdate;
  uint32_t                               lookaheadDepth = 1;
  std::vector<SwapOrMove>                swapCloseByLookahead;
  std::vector<std::pair<SwapOrMove, fp>> moveExactLookahead;
  MapperParameters                       parameters;
  std::deque<std::set<HwQubit>>          lastBlockedQubits;
  std::deque<AtomMove>                   lastMoves;
  Swap                                   lastSwap = Swap(0, 0);
  std::vector<fp>                        decayWeights;
  uint32_t                               nSwaps          = 0;
  uint32_t                               nMoves          = 0;
  fp                                     swapCloseByCost = 1;

  //  NeutralAtomMappingResults   results;
  HardwareQubits hardwareQubits;
  Mapping        mapping;
  bool           verbose = true;

  // Methods for layer creation
  void        createFrontLayer();
  void        updateFrontLayerByGate(GateList& gatesToExecute);
  void        updateFrontLayerByQubit();
  void        updateFrontLayerByCandidates();
  void        findFrontCandidates();
  void        updateLookaheadLayerByQubit();
  void        findLookaheadCandidates();
  void        updateLookaheadLayerByCandidates();
  void        mapGate(const Operation* op);
  static bool commutesWithAtQubit(const GateList&  layer,
                                  const Operation* opPointer,
                                  const Qubit&     qubit);
  static bool commuteAtQubit(const Operation* opPointer1,
                             const Operation* opPointer2, const Qubit& qubit);
  bool        isExecutable(const Operation* opPointer);
  void        addToFrontLayer(const Operation* opPointer);
  void        addToLookaheadLayer(const Operation* opPointer);
  void        reassignGatesToLayers();

  // Methods for estimation
  std::pair<uint32_t, fp> estimateNumSwapGates(const Operation* opPointer);
  std::pair<uint32_t, fp> estimateNumMove(const Operation* opPointer);
  bool                    swapGateBetter(const Operation* opPointer);

  // Methods for mapping
  Swap           findBestSwap();
  GateList       getExecutableGates();
  std::set<Swap> getAllPossibleSwaps();
  AtomMove       findBestAtomMove(const bool gateExecuted);
  MoveCombs      getMoveAwayCombinations(CoordIndex start, CoordIndex target);
  MultiQubitMovePos getMovePositionRec(MultiQubitMovePos   currentPos,
                                       const CoordIndices& gateCoords,
                                       const size_t&       maxNMoves);
  MoveCombs         getBestPossibleMoveCombinations(bool gateExecuted);
  MultiQubitMovePos getBestMovePos(const CoordIndices& gateCoords);

  void updateMapping(Swap swap);
  void updateMappingMove(AtomMove move);

  // Methods cost functions
  fp       distanceCost(const Swap& swap);
  void     initSwapAndMove(const GateList&                         layer,
                           std::vector<SwapOrMove>&                swapCloseBy,
                           std::vector<std::pair<SwapOrMove, fp>>& moveExact);
  fp       distancePerLayer(const Swap&                                   swap,
                            const std::vector<SwapOrMove>&                swapCloseBy,
                            const std::vector<std::pair<SwapOrMove, fp>>& moveExact);
  fp       moveCost(const AtomMove& move);
  fp       moveCostComb(const MoveComb& moveComb);
  fp       moveDistancePerLayer(const AtomMove& move, GateList& layer);
  fp       parallelMoveCost(const AtomMove& move);
  HwQubits getBestMultiQubitPosition(const Operation* opPointer);
  HwQubits getBestMultiQubitPositionRec(HwQubits remainingGateQubits,
                                        std::vector<HwQubit> selectedQubits,
                                        HwQubits remainingNearbyQubits);
  std::vector<std::pair<SwapOrMove, fp>>
            getExactMoveToPosition(const Operation* op, HwQubits position);
  MoveCombs getMoveCombinationsToPosition(HwQubits&                gateQubits,
                                          std::vector<CoordIndex>& position,
                                          bool gateExecuted);
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
  QuantumComputation map(qc::QuantumComputation& qc,
                         InitialMapping initialMapping, bool verbose = true);
  QuantumComputation mapAod(qc::QuantumComputation& qc);
  void setParameters(const MapperParameters& p) { this->parameters = p; }
};

} // namespace qc
