// Benchmark was created by MQT Bench on 2023-06-29
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: v1.0.0
// Qiskit version: {'qiskit-terra': '0.24.1', 'qiskit-aer': '0.12.0', 'qiskit-ignis': None, 'qiskit-ibmq-provider': '0.20.2', 'qiskit': '0.43.1', 'qiskit-nature': '0.6.2', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.5.0', 'qiskit-machine-learning': '0.6.1'}
// Used Gate Set: ['rx', 'rz', 'cz', 'measure']

OPENQASM 2.0;
include "qelib1.inc";
qreg q[10];
creg meas[10];
rz(pi/2) q[0];
rx(pi/2) q[0];
rz(pi/2) q[1];
rx(pi/2) q[1];
rz(pi/2) q[2];
rx(pi/2) q[2];
rz(pi/2) q[3];
rx(pi/2) q[3];
rz(pi/2) q[4];
rx(pi/2) q[4];
rz(pi/2) q[5];
rx(pi/2) q[5];
rz(pi/2) q[6];
rx(pi/2) q[6];
rz(pi/2) q[7];
rx(pi/2) q[7];
rz(pi/2) q[8];
rx(pi/2) q[8];
rz(-pi/2) q[9];
rx(-pi/2) q[9];
rz(-pi/2) q[9];
cz q[9],q[8];
rx(-pi/2) q[8];
rz(-pi/2) q[8];
cz q[8],q[7];
rx(-pi/2) q[7];
rz(-pi/2) q[7];
cz q[7],q[6];
rx(-pi/2) q[6];
rz(-pi/2) q[6];
cz q[6],q[5];
rx(-pi/2) q[5];
rz(-pi/2) q[5];
cz q[5],q[4];
rx(-pi/2) q[4];
rz(-pi/2) q[4];
cz q[4],q[3];
rx(-pi/2) q[3];
rz(-pi/2) q[3];
cz q[3],q[2];
rx(-pi/2) q[2];
rz(-pi/2) q[2];
cz q[2],q[1];
rx(-pi/2) q[1];
rz(-pi/2) q[1];
cz q[1],q[0];
rx(-pi/2) q[0];
rz(-pi/2) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
measure q[4] -> meas[4];
measure q[5] -> meas[5];
measure q[6] -> meas[6];
measure q[7] -> meas[7];
measure q[8] -> meas[8];
measure q[9] -> meas[9];