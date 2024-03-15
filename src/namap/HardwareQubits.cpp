//
// Created by Ludwig Schmid on 19.10.23.
//

#include "namap/HardwareQubits.hpp"

#include "namap/Mapping.hpp"

namespace qc {

void HardwareQubits::initializeFromGraph(const NeutralAtomArchitecture& arch, const DAG& dag){
  auto archNqubits    = arch.getNqubits();
  auto archNpositions = arch.getNpositions();
  auto archNrows      = arch.getNrows();
  auto archNcolumns   = arch.getNcolumns();

  std::vector<CoordIndex> qubitIndices(archNqubits, std::numeric_limits<unsigned int>::max());
  std::vector<CoordIndex> hwIndices(archNpositions, std::numeric_limits<unsigned int>::max());

  // interaction graph
  std::vector<std::vector<double>> circGraph(dag.size(), std::vector<double> (dag.size(), 0.0) );
  std::vector<std::pair<int, std::pair<int, double>>> circGraph_degree(dag.size());
  std::vector<std::vector<std::pair<int,double>>> circGraph_neighbor(dag.size());
  for(uint32_t qubit=0; qubit<dag.size(); ++qubit){
    for(auto opPtr : dag[qubit]){
        auto* op = opPtr->get();
        if(op->getUsedQubits().size() > 1){
            for(auto i : op->getUsedQubits()){
                if(i!=qubit){
                    circGraph[qubit][i] += 1;
                }
            }
        }
    }
  }
  // print interaction graph
  /*
  std::cout << std::endl << " # circGraph table " << std::endl;
  std::cout << "i↓  j→";
  for(uint32_t i=0; i<dag.size(); i++)
      std::cout << std::setw(4) << i;
  std::cout << std::endl;
  for(uint32_t i=0; i<dag.size(); i++)
  {
      std::cout << std::setw(6) << i;
      for(uint32_t j=0; j<dag.size(); j++)
          std::cout << std::setw(4) << circGraph[i][j];
      std::cout << std::endl;
  }
  std::cout << std::endl;           
  */
  // generate graph matching queue
  for(uint32_t qubit=0; qubit<dag.size(); qubit++){
    int cnt = 0;
    double sum = 0;
    for(uint32_t i=0; i<dag.size(); i++){
        double weight = circGraph[qubit][i];
        if(weight > 0){
            cnt++;
            sum += weight;
            circGraph_neighbor[qubit].emplace_back( i, weight );
        }
    }
    circGraph_degree[qubit] = std::make_pair( qubit, std::make_pair(cnt, sum));
  }
  sort(circGraph_degree.begin(), circGraph_degree.end(), 
        [](std::pair<int, std::pair<int,double>> a, std::pair<int, std::pair<int, double>> b){
            if(a.second.first == b.second.first)
                return a.second.second > b.second.second;
            else
                return a.second.first > b.second.first;
        });
  std::queue<uint32_t> circGraph_queue;
  for(const auto i : circGraph_degree){
    circGraph_queue.push(i.first);
  }
  for(auto& innerVec : circGraph_neighbor){
    sort(innerVec.begin(), innerVec.end(), 
        [](std::pair<int,double>& a, std::pair<int,double>& b){
            return a.second > b.second;
        });
  }
  // graph matching
  //std::cout << "## graph_matching: " << std::endl;
  bool firstCenter = true;
  uint32_t nMapped = 0;
  uint32_t archCenterX = (archNcolumns % 2 == 0) ? (archNcolumns/2 - 1) : (archNcolumns - 1)/2;
  uint32_t archCenterY = (archNrows % 2 == 0) ? (archNrows/2 - 1) : (archNrows - 1)/2;
  uint32_t archCenter = archCenterY * archNcolumns + archCenterX;
  //std::cout << "#HwQubit: " << archNqubits << "-> nRows: " << archNrows << ", nCols: " << archNcolumns << std::endl;
  //std::cout << "center: " << archCenter << " = (" << archCenterX << ", " << archCenterY << ")" << std::endl;

  while(!circGraph_queue.empty() && nMapped!=dag.size()){
    uint32_t qc = circGraph_queue.front();
    uint32_t hc = std::numeric_limits<unsigned int>::max(); //hardwarCenter
    //std::cout << "* qc: " << qc << std::endl;
    // center mappingx
    if(firstCenter){
        hc = archCenter;
        qubitIndices[qc] = hc;
        hwIndices[hc] = qc;
        firstCenter = false;
        nMapped++;
    }
    else if(qubitIndices[qc]==std::numeric_limits<unsigned int>::max()){
        //find position
        //std::cout << "need position for qc!" << std::endl;

        //ref loc
        std::vector<int> refLoc;
        for(auto i : circGraph_neighbor[qc]){
            if(qubitIndices[i.first] != std::numeric_limits<unsigned int>::max()){
                refLoc.push_back( qubitIndices[i.first] );
            }
        }
        /*
        std::cout << "refLoc: ";
        for(auto i : refLoc){
            std::cout << i << ", ";
        }
        std::cout << std::endl;
        */
        //candidate loc
        std::vector< std::pair<int,int> > distCandiLoc;
        std::vector<int> initCandiLoc;
        for(int i=0; i < archNpositions; i++){
            if(hwIndices[i] == std::numeric_limits<unsigned int>::max()){
                initCandiLoc.push_back(i);
            }
        }
        for(auto v : initCandiLoc){
            int distSum = 0;
            for(auto r : refLoc){
                int dist = std::abs(v%archNcolumns - r%archNcolumns) + std::abs(v/archNcolumns-r/archNcolumns);
                distSum += dist;
            }
            distCandiLoc.emplace_back(v, distSum);
        }
        sort(distCandiLoc.begin(), distCandiLoc.end(), 
            [](std::pair<int,int> a, std::pair<int,int> b){
                return a.second < b.second;
        });
        /*
        std::cout << "candiLoc: ";
        for(auto i : distCandiLoc){
            std::cout << i.first << ", ";
        }
        std::cout << std::endl;
        */
        hc = distCandiLoc[0].first;
        //std::cout << "hc: " << hc << std::endl;
        qubitIndices[qc] = hc;
        hwIndices[hc] = qc;
        nMapped++;
    }
    else{
        hc = qubitIndices[qc];
    }

    // neighbor mapping
    if(!circGraph_neighbor[qc].empty()){
        int idx_qc_n = 0;
        std::vector<int> qc_n;
        //std::cout << "q_neighbor: ";
        for(auto i : circGraph_neighbor[qc]){
            if( qubitIndices[i.first] != std::numeric_limits<unsigned int>::max() ) continue;
            if( idx_qc_n >= 4 ) continue;
            else{
                qc_n.push_back( i.first );
                idx_qc_n++;
                //std::cout << i.first << ", ";
            }
        }
        //std::cout << std::endl;

        std::vector<int> hw_n;
        if(hc!=std::numeric_limits<unsigned int>::max() && qc_n.size()>0){
            if((hc+1)%archNcolumns != 0 && hwIndices[hc+1]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc+1); //right
            if(hc/archNcolumns < (archNrows-1) && hwIndices[hc+archNcolumns]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc+archNcolumns);  //down
            if(hc%archNcolumns != 0 && hwIndices[hc-1]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc-1); //left
            if(hc>archNcolumns && hwIndices[hc-archNcolumns]==std::numeric_limits<unsigned int>::max()) hw_n.push_back(hc-archNcolumns); //up
            /*
            std::cout << "hw_neighbor: ";
            for(auto i : hw_n){
                std::cout << i << ", ";
            }
            std::cout << std::endl;
            */
        }

        int minSize = std::min( qc_n.size(), hw_n.size() );
        for(int i=0; i<minSize; i++){
            int qc_i = qc_n[i];
            int hw_i = hw_n[i];
            qubitIndices[qc_i] = hw_i;
            hwIndices[hw_i] = qc_i;
            nMapped++;
        }

    }
    circGraph_queue.pop();
  }
  
  int hwIndex = 0;
  for (uint32_t i = 0; i < archNqubits; ++i) {
    if(qubitIndices[i] == std::numeric_limits<unsigned int>::max()){
        if(hwIndices[hwIndex] != std::numeric_limits<unsigned int>::max()){
            do{
                hwIndex += 1;
            }while(hwIndices[hwIndex] != std::numeric_limits<unsigned int>::max());
        }
        hwToCoordIdx.insert({i, hwIndex});
        //std::cout << "hwToCoordIdx{" << i << ", " << hwIndex << ")" << std::endl;
        hwIndex++;
    }
    else{
        hwToCoordIdx.insert({i, qubitIndices[i]});
        //std::cout << "hwToCoordIdx{" << i << ", " << qubitIndices[i] << ")" << std::endl;
    }
  }
  
  swapDistances = SymmetricMatrix(archNqubits, -1);
}

void HardwareQubits::initTrivialSwapDistances() {
  swapDistances = SymmetricMatrix(arch.getNqubits());
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      swapDistances(i, j) =
          arch.getSwapDistance(hwToCoordIdx.at(i), hwToCoordIdx.at(j));
    }
  }
}

void HardwareQubits::initNearbyQubits() {
  for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
    computeNearbyQubits(i);
  }
}

void HardwareQubits::computeSwapDistance(HwQubit q1, HwQubit q2) {
  std::queue<HwQubit>  q;
  std::vector<bool>    visited(swapDistances.getSize(), false);
  std::vector<HwQubit> parent(swapDistances.getSize(), q2);

  q.push(q1);
  visited[q1] = true;
  parent[q1]  = q1;
  bool found  = false;
  while (!q.empty() && !found) {
    auto current = q.front();
    q.pop();
    for (const auto& nearbyQubit : nearbyQubits.at(current)) {
      if (!visited[nearbyQubit]) {
        q.push(nearbyQubit);
        visited[nearbyQubit] = true;
        parent[nearbyQubit]  = current;
        if (nearbyQubit == q2) {
          found = true;
          break;
        }
      }
    }
  }
  if (!found) {
    swapDistances(q1, q2) = std::numeric_limits<fp>::infinity();
    return;
  }
  // recreate path
  std::vector<HwQubit> path;
  auto                 current = q2;
  while (current != q1) {
    path.emplace_back(current);
    current = parent[current];
  }
  path.emplace_back(q1);
  // update swap distances along path
  for (uint32_t start = 0; start < path.size() - 1; ++start) {
    for (uint32_t end = start + 1; end < path.size(); ++end) {
      swapDistances(path[start], path[end]) = end - start - 1;
    }
  }
  //
  //  swapDistances(q1, q2) = static_cast<fp>(path.size() - 2);
}

void HardwareQubits::resetSwapDistances() {
  // TODO Improve to only reset the swap distances necessary (use a breadth
  // first search)
  swapDistances = SymmetricMatrix(arch.getNqubits(), -1);
}

void HardwareQubits::move(HwQubit hwQubit, CoordIndex newCoord) {
  if (newCoord >= arch.getNpositions()) {
    throw std::runtime_error("Invalid coordinate");
  }
  // check if new coordinate is already occupied
  for (const auto& [qubit, coord] : hwToCoordIdx) {
    if (coord == newCoord) {
      throw std::runtime_error("Coordinate already occupied");
    }
  }

  // remove qubit from old nearby qubits
  auto prevNearbyQubits = nearbyQubits.at(hwQubit);
  for (const auto& qubit : prevNearbyQubits) {
    nearbyQubits.at(qubit).erase(std::find(
        nearbyQubits.at(qubit).begin(), nearbyQubits.at(qubit).end(), hwQubit));
  }
  // move qubit and compute new nearby qubits
  hwToCoordIdx.at(hwQubit) = newCoord;
  computeNearbyQubits(hwQubit);

  // add qubit to new nearby qubits
  auto newNearbyQubits = nearbyQubits.at(hwQubit);
  for (const auto& qubit : newNearbyQubits) {
    nearbyQubits.at(qubit).insert(hwQubit);
  }

  // update/reset swap distances
  resetSwapDistances();
}

std::vector<Swap> HardwareQubits::getNearbySwaps(qc::HwQubit q) const {
  std::vector<Swap> swaps;
  swaps.reserve(nearbyQubits.size());
  for (const auto& nearbyQubit : nearbyQubits.at(q)) {
    swaps.emplace_back(q, nearbyQubit);
  }
  return swaps;
}

void HardwareQubits::computeNearbyQubits(qc::HwQubit q) {
  std::set<HwQubit> newNearbyQubits;
  auto              coordQ = hwToCoordIdx.at(q);
  for (const auto& coord : hwToCoordIdx) {
    if (coord.first == q) {
      continue;
    }
    if (arch.getEuclidianDistance(coordQ, coord.second) <=
        arch.getInteractionRadius()) {
      newNearbyQubits.insert(coord.first);
    }
  }
  nearbyQubits.insert_or_assign(q, newNearbyQubits);
}

fp HardwareQubits::getAllToAllSwapDistance(std::set<HwQubit>& qubits) {
  // two qubit gates
  if (qubits.size() == 2) {
    auto it = qubits.begin();
    auto q1 = *it;
    auto q2 = *(++it);
    return getSwapDistance(q1, q2);
  }
  // for n > 2 all qubits need to be within the interaction radius of each other
  fp totalDistance = 0;
  for (auto it1 = qubits.begin(); it1 != qubits.end(); ++it1) {
    for (auto it2 = std::next(it1); it2 != qubits.end(); ++it2) {
      totalDistance += getSwapDistance(*it1, *it2);
    }
  }
  return totalDistance;
}

std::set<HwQubit>
HardwareQubits::getBlockedQubits(const std::set<HwQubit>& qubits) {
  std::set<HwQubit> blockedQubits;
  for (const auto& qubit : qubits) {
    for (uint32_t i = 0; i < arch.getNqubits(); ++i) {
      if (i == qubit) {
        continue;
      }
      // TODO improve by using the nearby coords as a preselection
      auto const distance =
          arch.getEuclidianDistance(hwToCoordIdx.at(qubit), hwToCoordIdx.at(i));
      if (distance <= arch.getBlockingFactor() * arch.getInteractionRadius()) {
        blockedQubits.insert(i);
      }
    }
  }
  return blockedQubits;
}

std::set<CoordIndex> HardwareQubits::getNearbyFreeCoordinates(HwQubit q) {
  return getNearbyFreeCoordinatesByCoord(hwToCoordIdx.at(q));
}

std::set<CoordIndex>
HardwareQubits::getNearbyFreeCoordinatesByCoord(qc::CoordIndex idx) {
  std::set<CoordIndex> nearbyFreeCoordinates;
  for (auto const& coordIndex : this->arch.getNearbyCoordinates(idx)) {
    if (!this->isMapped(coordIndex)) {
      nearbyFreeCoordinates.insert(coordIndex);
    }
  }
  return nearbyFreeCoordinates;
}

std::set<CoordIndex>
HardwareQubits::getNearbyOccupiedCoordinates(qc::HwQubit q) {
  auto nearbyHwQubits = this->getNearbyQubits(q);
  return this->getCoordIndices(nearbyHwQubits);
}

std::set<CoordIndex>
HardwareQubits::getNearbyOccupiedCoordinatesByCoord(qc::CoordIndex idx) {
  auto nearbyHwQubits = this->getNearbyQubits(this->getHwQubit(idx));
  return this->getCoordIndices(nearbyHwQubits);
}

std::vector<CoordIndex>
HardwareQubits::findClosestFreeCoord(CoordIndex coord, Direction direction,
                                     const CoordIndices& excludeCoord) {
  // return the closest free coord in general
  // and the closest free coord in the given direction
  std::vector<CoordIndex> closestFreeCoords;
  std::queue<CoordIndex>  queue;
  queue.push(coord);
  std::set<CoordIndex> visited;
  visited.insert(coord);
  bool foundClosest = false;
  while (!queue.empty()) {
    auto currentCoord = queue.front();
    queue.pop();
    auto nearbyCoords = this->arch.getNN(currentCoord);
    for (const auto& nearbyCoord : nearbyCoords) {
      if (std::find(visited.rbegin(), visited.rend(), nearbyCoord) ==
          visited.rend()) {
        visited.insert(nearbyCoord);
        if (!this->isMapped(nearbyCoord) &&
            std::find(excludeCoord.begin(), excludeCoord.end(), nearbyCoord) ==
                excludeCoord.end()) {
          if (!foundClosest) {
            closestFreeCoords.push_back(nearbyCoord);
          }
          foundClosest = true;
          if (direction == arch.getVector(coord, nearbyCoord).direction) {
            closestFreeCoords.push_back(nearbyCoord);
            return closestFreeCoords;
          }
        } else {
          queue.push(nearbyCoord);
        }
      }
    }
  }
  return closestFreeCoords;
}

[[maybe_unused]] fp HardwareQubits::getSwapDistanceMove(CoordIndex idx,
                                                        HwQubit    target) {
  auto nearbyCoords    = this->arch.getNearbyCoordinates(idx);
  fp   minSwapDistance = std::numeric_limits<fp>::infinity();
  if (nearbyCoords.find(this->getCoordIndex(target)) != nearbyCoords.end()) {
    return 0;
  }
  for (const auto& nearbyCoord : nearbyCoords) {
    if (this->isMapped(nearbyCoord)) {
      auto swapDistance =
          this->getSwapDistance(target, this->getHwQubit(nearbyCoord)) + 1;
      if (swapDistance < minSwapDistance) {
        minSwapDistance = swapDistance;
      }
    }
  }
  return minSwapDistance;
}

} // namespace qc
