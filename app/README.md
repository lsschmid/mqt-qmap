# COMMENTS

This is the source for a simple `main.cpp` to use the Hybrid Mapper proposed in the paper.
In the following some more information on how to use the mapper to reproduce the results of the paper und the corresponding experiments is given.

## WIP !!!

The mapper is still work in progress and not yet ready for general use.
An updated/rewritten version of the mapper will be made available in the future as part of the MQT family.

## How to run the Hybrid Mapper

1. Put the input quantum circuits (.qasm2) in the `input` folder.
2. Run the app with `NaMain runIdx input_directory output_directory lookaheadGate lookaheadShuttling gateDecay shuttlingTimeWeight gateWeight shuttlingWeight verbose json_config_file_path initialCoordinateMapping initialCircuitMapping`
   - `runIdx` is an integer that is used to identify the run. It is used to create the output folder.
   - `input_directory` is the path to the input directory.
   - `output_directory` is the path to the output directory.
   - `lookaheadGate` is the lookahead weight for the gate-based mapping: $w_{l,\mathrm{g}}$.
   - `lookaheadShuttling` is the lookahead for the shuttling-based mapping: $w_{l,\mathrm{s}}$.
   - `gateDecay` is the decay factor for the gate placement algorithm $\lambda^t$.
   - `shuttlingTimeWeight` is the weight for the time weight in the shuttling objective function $w^t$.
   - `gateWeight` is the weight for the gate-based mapping: $\alpha_\mathrm{g}$.
   - `shuttlingWeight` is the weight for the shuttling-based mapping: $\alpha_\mathrm{s}$.
   - `verbose` is a boolean that indicates whether the app should print additional information.
   - `json_config_file_path` is the path to the JSON configuration file that contains the device information.
   - `initialCoordinateMapping` is an identifier for the atom mapping $f_\mathrm{a}$ (only "trivial/identity" supported at the moment).
   - `initialCircuitMapping` is an identifier for the qubit mapping $f_\mathrm{q}$ (only "trivial/identity" supported at the moment).
3. The output is written to the `output` folder. The output folder is created if it does not exist. It contains the following files
   - parameters_idxRun.csv : file saving the run parameters
   - "circuit.qasm_ext" : file saving the circuit including the SWAP and shuttling MOVE operations
   - "circuit.qasm_aod" : file saving the circuit with SWAPs decomposed to CZ + H gates and shuttling MOVE operations scheduled and converted to AOD operations (activate, deactivate, move).
   - idxRun.csv : file saving the run results (mapping fidelity and qubit idle times).

For further information, please check out directly the code of `main.cpp` and the corresponding libraries.
For any questions regarding the code please contact the authors of the paper.

## Data and evaluation results

The folder `numeric_evaluations` contains the data and the scripts used for the numeric evaluations presented in the paper.
It contains the following subfolders:

- `benchmarks` : contains the benchmark circuits used for the numeric evaluations.
- `json` : contains the JSON configuration files for the numeric evaluations (device information).
- `output_evaluations` : contains the output (descirbed above) of the numeric evaluations for the different hardware configurations and the different mapper configurations (named: `gateWeight_shuttlingWeight`).
