# paramters: runIdx inputDir outputDir lookaheadGate lookaheadShuttling decay shuttlingTimeWeight shuttlingMakeExecutalbeBonus shuttlingMultiQubitWeight GateMultiQubitWeight gateWeight ShuttlingWeight verbose jsonFile initCoordMapping initCircuitMapping
#./NaMain 1 benchmarks output 0.1 0.1 0 0.1 1 1 1 1 1 0 rubidium.json trivial identity

#./NaMain 1 benchmarks_multi_qubit output 0.1 0.1 0 0.1 1 1 1 1 1 1 rubidium.json trivial identity
./NaMain 1 benchmarks/50-qubit output 0.1 0.1 0 0.1 1 1 ./json/rubidium_shuttling.json trivial identity
