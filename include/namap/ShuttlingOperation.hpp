//
// Created by Ludwig Schmid on 05.10.23.
//
#pragma once

#include "operations/Operation.hpp"

namespace qc {

enum ShuttlingType : std::uint8_t {
  Move,
  Activate,
  Deactivate,
};

const inline static std::unordered_map<std::string, qc::ShuttlingType>
    SHUTTLING_OP_NAME_TO_TYPE = {
        {"move", ShuttlingType::Move},
        {"activate", ShuttlingType::Activate},
        {"deactivate", ShuttlingType::Deactivate},
};

class ShuttlingOperation : public Operation {};

class AODOperation : public Operation {};

} // namespace qc
