//
// Created by Ludwig Schmid on 19.10.23.
//

#include "namap/Mapping.hpp"

namespace qc {
void Mapping::swap(Swap swap) {
  auto q1 = swap.first;
  auto q2 = swap.second;
  if (this->isMapped(q1) && this->isMapped(q2)) {
    auto circQ1 = this->getCircQubit(q1);
    auto circQ2 = this->getCircQubit(q2);
    this->removeCircuitQubit(circQ1);
    this->removeCircuitQubit(circQ2);
    this->setCircuitQubit(circQ2, q1);
    this->setCircuitQubit(circQ1, q2);
  } else if (this->isMapped(q1) && !this->isMapped(q2)) {
    this->setCircuitQubit(this->getCircQubit(q1), q2);
    this->removeCircuitQubit(q1);
  } else if (this->isMapped(q2) && !this->isMapped(q1)) {
    this->setCircuitQubit(this->getCircQubit(q2), q1);
    this->removeCircuitQubit(q2);
  } else {
    throw std::runtime_error("Cannot swap unmapped qubits");
  }
}

} // namespace qc
