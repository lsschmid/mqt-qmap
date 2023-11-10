//
// Created by Ludwig Schmid on 19.10.23.
//

#pragma once

#include "operations/AodOperation.hpp"
#include "utils.hpp"

namespace qc {
using HwQubit      = uint32_t;
using Swap         = std::pair<HwQubit, HwQubit>;
using CoordIndex   = std::uint32_t;
using CoordIndices = std::vector<CoordIndex>;
using AtomMove     = std::pair<CoordIndex, CoordIndex>;
using Qubits       = std::set<Qubit>;
using HwQubits     = std::set<HwQubit>;
using HwPositions  = std::vector<HwQubits>;
using SwapOrMove   = std::pair<HwQubit, HwQubit>;
using OpPointer    = std::unique_ptr<AodOperation>;

} // namespace qc
