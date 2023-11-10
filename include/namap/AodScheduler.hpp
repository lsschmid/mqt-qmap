//
// Created by Ludwig Schmid on 30.10.23.
//

#pragma once

#include "QuantumComputation.hpp"
#include "namap/NeutralAtomArchitecture.hpp"

#include <utility>

namespace qc {
enum class ActivationMerge { Impossible, Trivial, Merge, Append };

class AodScheduler {
protected:
  struct AodActivationHelper {
    struct AodMove {
      uint32_t init;
      fp       delta;
      int32_t  offset;

      AodMove(uint32_t init, fp delta, int32_t offset)
          : init(init), delta(delta), offset(offset) {}
    };
    struct AodActivation {
      // first: x, second: delta x, third: offset x
      std::vector<AodMove*> activateXs;
      std::vector<AodMove*> activateYs;
      std::vector<AtomMove> moves;

      AodActivation(const AodMove& activateXs, const AodMove& activateYs,
                    const AtomMove& move)
          : activateXs({new AodMove(activateXs)}),
            activateYs({new AodMove(activateYs)}), moves({move}) {}
      AodActivation(const Dimension dim, const AodMove& activate,
                    const AtomMove& move)
          : moves({move}) {
        if (dim == Dimension::X) {
          activateXs = {new AodMove(activate)};
        } else {
          activateYs = {new AodMove(activate)};
        }
      }

      [[nodiscard]] std::vector<AodMove*> inline getActivates(
          Dimension dim) const {
        if (dim == Dimension::X) {
          return activateXs;
        }
        return activateYs;
      }
    };

    NeutralAtomArchitecture    arch;
    std::vector<AodActivation> allActivations;
    OpType                     type;

    [[nodiscard]] std::vector<AodMove*> getAodMoves(Dimension dim,
                                                    uint32_t  x) const;

    std::pair<ActivationMerge, ActivationMerge>
    canAddActivation(const Coordinate& origin, const AtomMove& move,
                     MoveVector v) const;
    ActivationMerge canMergeActivation(Dimension dim, const Coordinate& origin,
                                       MoveVector v) const;
    void        addActivation(std::pair<ActivationMerge, ActivationMerge> merge,
                              const Coordinate& origin, const AtomMove& move,
                              MoveVector v);
    void        mergeActivation(Dimension dim, const AodActivation& activation);
    static void reAssignOffsets(std::vector<AodMove*>& aodMoves, int32_t sign);

    [[nodiscard]] uint32_t getMaxOffset(Dimension dim, uint32_t x,
                                        int32_t sign) const;

    [[nodiscard]] bool checkIntermediateSpace(Dimension dim, uint32_t x,
                                              int32_t sign) const;

    [[nodiscard]] static std::pair<AodOperation, AodOperation>
    getAodOperation(const AodActivation&           activation,
                    const NeutralAtomArchitecture& arch, OpType type);
    [[nodiscard]] std::vector<AodOperation> getAodOperations() const;

    AodActivationHelper(NeutralAtomArchitecture arch, OpType type)
        : arch(std::move(arch)), type(type) {}
  };

  struct MoveGroup {
    NeutralAtomArchitecture                    arch;
    std::vector<std::pair<AtomMove, uint32_t>> moves;
    std::vector<AodOperation>                  processedOpsInit;
    std::vector<AodOperation>                  processedOpsFinal;
    AodOperation                               processedOpShuttle;
    std::vector<CoordIndex>                    targetQubits;
    std::vector<CoordIndex>                    qubitsUsedByGates;

    bool                          canAdd(const AtomMove& move);
    void                          add(const AtomMove& move, const uint32_t idx);
    [[nodiscard]] inline uint32_t getFirstIdx() const {
      return moves.front().second;
    }
    [[nodiscard]] inline bool isFirstOpSet() const { return !moves.empty(); }
    static bool parallelCheck(const MoveVector& v1, const MoveVector& v2);

    static AodOperation
    connectAodOperations(const std::vector<AodOperation>& opsInit,
                         const std::vector<AodOperation>& opsFinal);

    MoveGroup(NeutralAtomArchitecture arch) : arch(std::move(arch)) {}
  };

  NeutralAtomArchitecture arch;
  QuantumComputation      qcScheduled;
  std::vector<MoveGroup>  moveGroups;

  void initMoveGroups(QuantumComputation& qc);
  void processMoveGroups();

public:
  explicit AodScheduler(const NeutralAtomArchitecture& arch)
      : arch(arch), qcScheduled(arch.getNpositions()) {}

  QuantumComputation        schedule(QuantumComputation& qc);
  [[nodiscard]] inline auto getNMoveGroups() const { return moveGroups.size(); }
};

} // namespace qc
