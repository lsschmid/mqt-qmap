//
// Created by Ludwig Schmid on 05.10.23.
//

#include "filesystem"
#include "namap/NeutralAtomArchitecture.hpp"
#include "namap/NeutralAtomMapper.hpp"

#include "gtest/gtest.h"

TEST(NeutralAtomTest, architcure) {
  std::filesystem::path const testDirPath = std::filesystem::current_path();
  std::string const           testDir     = testDirPath.string();
  auto filename                    = std::string(testDir) + "/rubidium.json";
  qc::NeutralAtomArchitecture arch = qc::NeutralAtomArchitecture(filename);
  EXPECT_EQ(arch.getCoordinate(0).getX(), 0);
  EXPECT_EQ(arch.getNcolumns(), 5);
  //  auto d = arch.getDistance(0, 1);
  auto d =
      arch.getEuclidianDistance(arch.getCoordinate(0), arch.getCoordinate(1));
  EXPECT_EQ(d, 1);
  EXPECT_EQ(arch.getSwapDistance(0, 1), 0);
  EXPECT_EQ(arch.getSwapDistance(0, 2), 0);
  EXPECT_EQ(arch.getSwapDistance(0, 3), 1);
}

TEST(NeutralAtomTest, mapper) {
  std::filesystem::path const testDirPath = std::filesystem::current_path();
  std::string const           testDir     = testDirPath.string();
  auto filename                    = std::string(testDir) + "/rubidium.json";
  qc::NeutralAtomArchitecture arch = qc::NeutralAtomArchitecture(filename);
  qc::NeutralAtomMapper       mapper =
      qc::NeutralAtomMapper(arch, static_cast<qc::InitialCoordinateMapping>(0));
  //  auto qc = qc::QuantumComputation("../../examples/3_17_13.qasm");
  auto qc = qc::QuantumComputation("../../examples/qft_10.qasm");

  mapper.map(qc);

  EXPECT_EQ(0, 0);
}
