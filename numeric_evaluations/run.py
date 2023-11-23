# run.py for neutral atom mapper
from __future__ import annotations

import subprocess as sp

# define parameter
binaryName = "./NaMain"
inputDir = "benchmarks"
outputDir = "output_circuit/mixedHW/0.99_1"
jsonFile = "json/rubidium_mixed.json"

initCoordMapping = "trivial"
initCircuitMapping = "identity"

run_idx = 0
time_limit = 10000


# for experiment
result = []


## function define ##
def ExecuteCommand(
    runIdx, lookaheadGate, lookaheadShuttling, decay, shuttlingTimeWeight, gateWeight, ShuttlingWeight, verbose
):
    exeStr = f"{binaryName} {runIdx} {inputDir} {outputDir} {lookaheadGate} {lookaheadShuttling} {decay} {shuttlingTimeWeight} {gateWeight} {ShuttlingWeight} {verbose} {jsonFile} {initCoordMapping} {initCircuitMapping}"
    print(exeStr)
    return_code, stdout, error = RunCommand(exeStr, time_limit)
    # print(stdout)
    if return_code == 0:
        print("Successfully run!")
        result_for_case = {}
        result_for_case["benchmark_name"] = None
        result_for_case["nSwaps"] = None
        result_for_case["nMoves"] = None
        #
        result_for_case["totalExecutionTimes"] = None
        result_for_case["totalIdleTime"] = None
        result_for_case["totalGateFidelities"] = None
        result_for_case["totalFidelities"] = None
        result_for_case["nCZs"] = None
        result_for_case["runtime"] = None
        #
        result_for_case["runIdx"] = runIdx
        result_for_case["gateWeight"] = gateWeight
        result_for_case["ShuttlingWeight"] = ShuttlingWeight

        lines = stdout.decode("utf-8").split("\n")
        benchmark_name = None
        for line in lines:
            if line.startswith(f"Mapping {inputDir}/"):
                benchmark_name = line.split(f"Mapping {inputDir}/")[1].split(".qasm")[0]
                result_for_case["benchmark_name"] = benchmark_name
            if "nSwaps" in line:
                nSwaps = int(line.split(":")[1].strip())
                result_for_case["nSwaps"] = nSwaps
            if "nMoves" in line:
                nMoves = int(line.split(":")[1].strip())
                result_for_case["nMoves"] = nMoves
            if "totalExecutionTimes" in line:
                totalExecutionTimes = float(line.split(":")[1].strip())
                result_for_case["totalExecutionTimes"] = totalExecutionTimes
            if "totalIdleTime" in line:
                totalIdleTime = float(line.split(":")[1].strip())
                result_for_case["totalIdleTime"] = totalIdleTime
            if "totalGateFidelities" in line:
                totalGateFidelities = float(line.split(":")[1].strip())
                result_for_case["totalGateFidelities"] = totalGateFidelities
            if "totalFidelities" in line:
                totalFidelities = float(line.split(":")[1].strip())
                result_for_case["totalFidelities"] = totalFidelities
            if "totalnCZs" in line:
                nCZs = int(line.split(":")[1].strip())
                result_for_case["nCZs"] = nCZs
            if "* runtime" in line:
                runtime = int(line.split(":")[1].strip())
                result_for_case["runtime"] = runtime
                if result_for_case:
                    result_for_case["runIdx"] = runIdx
                    result_for_case["gateWeight"] = gateWeight
                    result_for_case["ShuttlingWeight"] = ShuttlingWeight
                    result.append(result_for_case)
                    result_for_case = {}

        if result_for_case:
            result.append(result_for_case)
    else:
        print(f"Error running the command. Return code: {return_code}")


def RunCommand(exeStr, time_limit):
    try:
        proc = sp.Popen(exeStr, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = proc.communicate(timeout=time_limit)
        return proc.returncode, stdout, stderr
    except sp.TimeoutExpired:
        print("Time out!")
        proc.terminate()
        proc.kill()
        return None, None, None


## execute ##
grouped_results = {}
sorted_results = []

# ./NaMain 0 benchmarks output 0.1 0.1 0 0.1 1 1 0 rubidium.json trivial identity
for gateWeight in [0.99]:
    for shuttlingWeight in [1]:
        if gateWeight != 0 or shuttlingWeight != 0:
            ExecuteCommand(run_idx, 0.1, 0.1, 0, 0.1, gateWeight, shuttlingWeight, 0)
            run_idx += 1

## sort ##
for res in result:
    benchmark_name = res["benchmark_name"]
    if benchmark_name not in grouped_results:
        grouped_results[benchmark_name] = []
    grouped_results[benchmark_name].append(res)

sorted_reuslts = []
for benchmark_name, results in grouped_results.items():
    sorted_results.extend(sorted(results, key=lambda x: x["totalFidelities"], reverse=True))

## Print the results ##
print()
print()
print(
    "---------------------------------------------------------------------- Result -------------------------------------------------------------------------"
)
print(
    "|                                                    ||        parameter      ||       mappingResults     ||     schedulerResults     ||    runtime   |"
)
print(
    "|                     benchmark                      ||index|  gateW |  shutW ||  nCZs  | nSwaps | nMoves || ExecTime  |  Fidelities  ||  runtime[ms] |"
)
print(
    "-------------------------------------------------------------------------------------------------------------------------------------------------------"
)
for res in sorted_results:
    print(
        "|",
        f"{res['benchmark_name']:<50}",
        "||",
        #
        f"{res['runIdx']:<3}",
        "|",
        f"{res['gateWeight']:<6}",
        "|",
        f"{res['ShuttlingWeight']:<6}",
        "||",
        #
        f"{res['nCZs']:<6}",
        "|",
        f"{res['nSwaps']:<6}",
        "|",
        f"{res['nMoves']:<6}",
        "||",
        #
        f"{res['totalExecutionTimes']:<9}",
        "|",
        f"{res['totalFidelities']:<12}",
        "||",
        #
        f"{res['runtime']:<12}",
        "|",
    )
print(
    "-------------------------------------------------------------------------------------------------------------------------------------------------------"
)
