//
// Created by Ludwig Schmid on 06.11.23.
//

#include "filesystem"
#include "namap/NeutralAtomMapper.hpp"
#include "namap/NeutralAtomScheduler.hpp"

int main(int argc, char* argv[]) {
  if (argc != 14) {
    std::cerr
        << "Usage: " << argv[0]
        << " <runIdx> <input_directory> <output_directory> "
           "<lookaheadGate> <lookaheadShuttling> <gateDecay> "
           "<shuttlingTimeWeight> "
           "<gateWeight> <shuttlingWeight> <verbose> <json_config_file_path> "
           "<initialCoordinateMapping> <initialCircuitMapping>\n";
    return 1;
  }

  int         runIdx                   = std::atoi(argv[1]);
  std::string input_directory          = argv[2];
  std::string output_directory         = argv[3];
  double      lookaheadGate            = std::stod(argv[4]);
  double      lookaheadShuttling       = std::stod(argv[5]);
  double      gateDecay                = std::stod(argv[6]);
  double      shuttlingTimeWeight      = std::stod(argv[7]);
  double      gateWeight               = std::stod(argv[8]);
  double      shuttlingWeight          = std::stod(argv[9]);
  bool        verbose                  = std::atoi(argv[10]) != 0;
  std::string json_config_file_path    = argv[11];
  std::string initialCoordinateMapping = argv[12];
  std::string initialCircuitMapping    = argv[13];

  // Check if the output directory exists and create it if it doesn't.
  if (!std::filesystem::exists(output_directory)) {
    if (!std::filesystem::create_directory(output_directory)) {
      std::cerr << "Failed to create the output directory.\n";
      return 1;
    }
  }

  // init Mappings
  qc::InitialCoordinateMapping initialCoordinateMappingEnum;
  if (initialCoordinateMapping == "trivial") {
    initialCoordinateMappingEnum = qc::InitialCoordinateMapping::Trivial;
  } else if (initialCoordinateMapping == "random") {
    initialCoordinateMappingEnum = qc::InitialCoordinateMapping::Random;
  } else {
    std::cerr << "Unknown initial coordinate mapping: "
              << initialCoordinateMapping << "\n";
    return 1;
  }
  qc::InitialMapping initialCircuitMappingEnum;
  if (initialCircuitMapping == "identity") {
    initialCircuitMappingEnum = qc::InitialMapping::Identity;
  } else {
    std::cerr << "Unknown initial circuit mapping: " << initialCircuitMapping
              << "\n";
    return 1;
  }

  // read files
  std::vector<std::string> qasmFiles;
  for (const auto& entry :
       std::filesystem::directory_iterator(input_directory)) {
    if (entry.is_regular_file() && entry.path().extension() == ".qasm") {
      qasmFiles.push_back(entry.path().string());
    }
  }

  // output file
  std::ofstream ofsResults(output_directory + "/" + std::to_string(runIdx) +
                           ".csv");

  for (const auto& qasmFile : qasmFiles) {
    // create arch
    qc::NeutralAtomArchitecture arch =
        qc::NeutralAtomArchitecture(json_config_file_path);
    // start mapping
    auto startTime = std::chrono::high_resolution_clock::now();
    qc::NeutralAtomMapper mapper =
        qc::NeutralAtomMapper(arch, initialCoordinateMappingEnum);
    qc::MapperParameters mapperParameters;
    mapperParameters.lookaheadWeightSwaps = lookaheadGate;
    mapperParameters.lookaheadWeightMoves = lookaheadShuttling;
    mapperParameters.decay                = gateDecay;
    mapperParameters.shuttlingTimeWeight  = shuttlingTimeWeight;
    mapperParameters.gateWeight           = gateWeight;
    mapperParameters.shuttlingWeight      = shuttlingWeight;
    mapper.setParameters(mapperParameters);

    std::cout << "Mapping " << qasmFile << "\n";
    qc::QuantumComputation qc = qc::QuantumComputation(qasmFile);
    auto          qcMapped = mapper.map(qc, initialCircuitMappingEnum, verbose);
    std::ofstream ofs(output_directory + "/" +
                      std::filesystem::path(qasmFile).filename().string() +
                      "_" + std::to_string(runIdx) + ".qasm_ext");
    qcMapped.dumpOpenQASM(ofs);
    auto          qcAodMapped = mapper.convertToAod(qcMapped);
    std::ofstream ofs_aod(output_directory + "/" +
                          std::filesystem::path(qasmFile).filename().string() +
                          "_" + std::to_string(runIdx) + ".qasm_aod");
    qcAodMapped.dumpOpenQASM(ofs_aod);
    auto endTime = std::chrono::high_resolution_clock::now();
    auto timeTaken = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
    // do the scheduling
    qc::NeutralAtomScheduler scheduler = qc::NeutralAtomScheduler(arch);
    auto schedulerResults = scheduler.schedule(qcAodMapped, verbose);

    // dump the results
    ofsResults << std::filesystem::path(qasmFile).filename().string() + ", " +
                      schedulerResults.toCsv() + "\n";

    // dump the parameters
    std::ofstream ofs_params(output_directory + "/parameters_" +
                             std::to_string(runIdx) + ".txt");
    ofs_params << "lookaheadGate: " << lookaheadGate << "\n";
    ofs_params << "lookaheadShuttling: " << lookaheadShuttling << "\n";
    ofs_params << "gateDecay: " << gateDecay << "\n";
    ofs_params << "shuttlingTimeWeight: " << shuttlingTimeWeight << "\n";
    ofs_params << "gateWeight: " << gateWeight << "\n";
    ofs_params << "shuttlingWeight: " << shuttlingWeight << "\n";
    ofs_params << "verbose: " << verbose << "\n";
    ofs_params << "json_config_file_path: " << json_config_file_path << "\n";
    ofs_params << "initialCoordinateMapping: " << initialCoordinateMapping
               << "\n";
    ofs_params << "initialCircuitMapping: " << initialCircuitMapping << "\n";
    // close the file
    ofs_params.close();
    std::cout << "* runtime: " << timeTaken << std::endl;
  }

  return 0;
}
