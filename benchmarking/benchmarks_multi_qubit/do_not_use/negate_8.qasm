OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
x q[0];
x q[1];
x q[2];
x q[3];
x q[4];
x q[5];
x q[6];
x q[7];
cccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7];
ccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6];
cccccx q[0], q[1], q[2], q[3], q[4], q[5];
ccccx q[0], q[1], q[2], q[3], q[4];
cccx q[0], q[1], q[2], q[3];
ccx q[0], q[1], q[2];
cx q[0], q[1];
x q[0];
x q[0];
cx q[0], q[1];
ccx q[0], q[1], q[2];
cccx q[0], q[1], q[2], q[3];
ccccx q[0], q[1], q[2], q[3], q[4];
cccccx q[0], q[1], q[2], q[3], q[4], q[5];
ccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6];
cccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7];
cccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7];
ccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6];
cccccx q[0], q[1], q[2], q[3], q[4], q[5];
ccccx q[0], q[1], q[2], q[3], q[4];
cccx q[0], q[1], q[2], q[3];
ccx q[0], q[1], q[2];
cx q[0], q[1];
x q[0];
x q[0];
cx q[0], q[1];
ccx q[0], q[1], q[2];
cccx q[0], q[1], q[2], q[3];
ccccx q[0], q[1], q[2], q[3], q[4];
cccccx q[0], q[1], q[2], q[3], q[4], q[5];
ccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6];
cccccccx q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7];
x q[0];
x q[1];
x q[2];
x q[3];
x q[4];
x q[5];
x q[6];
x q[7];