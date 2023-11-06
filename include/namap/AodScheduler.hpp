//
// Created by Ludwig Schmid on 30.10.23.
//

#pragma once

#include "QuantumComputation.hpp"
#include "namap/NeutralAtomArchitecture.hpp"

#include <utility>

namespace qc {
class AodScheduler {
protected:
  struct AodActivationHelper {
    struct AodMove {
      uint32_t init;
      fp       delta;
      int32_t  offset;
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

      [[nodiscard]] std::vector<AodMove*> inline getActivates(
          Dimension dim) const {
        if (dim == Dimension::X) {
          return activateXs;
        }
        return activateYs;
      }
      void addAodMove(Dimension dim, const AodMove& aodMove) {
        if (dim == Dimension::X) {
          activateXs.push_back(new AodMove(aodMove));
        } else {
          activateYs.push_back(new AodMove(aodMove));
        }
      }
    };

    NeutralAtomArchitecture    arch;
    std::vector<AodActivation> allActivations;
    uint32_t                   nIntermediateLevels;

    [[nodiscard]] std::vector<AodMove*> getAodMoves(Dimension dim,
                                                    uint32_t  x) const;

    bool        addActivation(const Coordinate& origin, const AtomMove& move,
                              MoveVector v);
    void        mergeActivation(Dimension dim, const AodActivation& activation);
    static void reAssignOffsets(std::vector<AodMove*>& aodMoves, int32_t sign);

    [[nodiscard]] uint32_t getMaxOffset(Dimension dim, uint32_t x,
                                        int32_t sign) const;

    [[nodiscard]] bool checkIntermediateSpace(Dimension dim, uint32_t x,
                                              int32_t sign) const;

    [[nodiscard]] static std::pair<OpPointer, OpPointer>
    getAodOperation(const AodActivation&           activation,
                    const NeutralAtomArchitecture& arch);
    [[nodiscard]] std::vector<OpPointer> getAodOperations() const;

    AodActivationHelper(NeutralAtomArchitecture arch) : arch(std::move(arch)) {}
  };

  struct MoveGroup {
    NeutralAtomArchitecture                    arch;
    std::vector<std::pair<AtomMove, uint32_t>> moves;
    std::vector<OpPointer>                     processedOps;

    bool                          canAdd(const AtomMove& move);
    void                          add(const AtomMove& move, const uint32_t idx);
    [[nodiscard]] inline uint32_t getFirstIdx() const {
      return moves.front().second;
    }
    [[nodiscard]] inline bool isFirstOpSet() const { return !moves.empty(); }
    static bool parallelCheck(const MoveVector& v1, const MoveVector& v2);

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

  QuantumComputation schedule(QuantumComputation& qc);
};

} // namespace qc
