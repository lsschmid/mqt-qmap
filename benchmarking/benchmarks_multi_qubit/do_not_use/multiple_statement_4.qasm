OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
cx q[5], q[9];
cx q[6], q[10];
cx q[7], q[11];
cx q[6], q[7];
cx q[5], q[6];
ccx q[4], q[8], q[5];
ccx q[5], q[9], q[6];
ccx q[6], q[10], q[7];
cx q[7], q[11];
ccx q[6], q[10], q[7];
cx q[6], q[10];
ccx q[5], q[9], q[6];
cx q[5], q[9];
ccx q[4], q[8], q[5];
cx q[4], q[8];
cx q[5], q[6];
cx q[6], q[7];
cx q[5], q[9];
cx q[6], q[10];
cx q[7], q[11];
cx q[8], q[12];
cx q[9], q[13];
cx q[10], q[14];
cx q[11], q[15];
x q[8];
x q[9];
x q[10];
x q[11];
cx q[5], q[9];
cx q[6], q[10];
cx q[7], q[11];
cx q[6], q[7];
cx q[5], q[6];
ccx q[4], q[8], q[5];
ccx q[5], q[9], q[6];
ccx q[6], q[10], q[7];
cx q[7], q[11];
ccx q[6], q[10], q[7];
cx q[6], q[10];
ccx q[5], q[9], q[6];
cx q[5], q[9];
ccx q[4], q[8], q[5];
cx q[4], q[8];
cx q[5], q[6];
cx q[6], q[7];
cx q[5], q[9];
cx q[6], q[10];
cx q[7], q[11];
x q[8];
x q[9];
x q[10];
x q[11];
x q[4];
x q[5];
x q[6];
x q[7];
cx q[5], q[13];
cx q[6], q[14];
cx q[7], q[15];
cx q[6], q[7];
cx q[5], q[6];
ccx q[4], q[12], q[5];
ccx q[5], q[13], q[6];
ccx q[6], q[14], q[7];
cx q[7], q[15];
ccx q[6], q[14], q[7];
cx q[6], q[14];
ccx q[5], q[13], q[6];
cx q[5], q[13];
ccx q[4], q[12], q[5];
cx q[4], q[12];
cx q[5], q[6];
cx q[6], q[7];
cx q[5], q[13];
cx q[6], q[14];
cx q[7], q[15];
x q[4];
x q[5];
x q[6];
x q[7];
x q[12];
x q[13];
x q[14];
x q[15];
cx q[13], q[9];
cx q[14], q[10];
cx q[15], q[11];
cx q[14], q[15];
cx q[13], q[14];
ccx q[8], q[12], q[13];
ccx q[9], q[13], q[14];
ccx q[10], q[14], q[15];
cx q[15], q[11];
ccx q[10], q[14], q[15];
cx q[14], q[10];
ccx q[9], q[13], q[14];
cx q[13], q[9];
ccx q[8], q[12], q[13];
cx q[12], q[8];
cx q[13], q[14];
cx q[14], q[15];
cx q[13], q[9];
cx q[14], q[10];
cx q[15], q[11];
x q[4];
x q[5];
x q[6];
x q[7];
cx q[5], q[13];
cx q[6], q[14];
cx q[7], q[15];
cx q[6], q[7];
cx q[5], q[6];
ccx q[4], q[12], q[5];
ccx q[5], q[13], q[6];
ccx q[6], q[14], q[7];
cx q[7], q[15];
ccx q[6], q[14], q[7];
cx q[6], q[14];
ccx q[5], q[13], q[6];
cx q[5], q[13];
ccx q[4], q[12], q[5];
cx q[4], q[12];
cx q[5], q[6];
cx q[6], q[7];
cx q[5], q[13];
cx q[6], q[14];
cx q[7], q[15];
x q[4];
x q[5];
x q[6];
x q[7];
x q[12];
x q[13];
x q[14];
x q[15];
cx q[12], q[8];
cx q[13], q[9];
cx q[14], q[10];
cx q[15], q[11];
x q[0];
x q[1];
x q[2];
x q[3];
cx q[9], q[1];
cx q[10], q[2];
cx q[11], q[3];
cx q[10], q[11];
cx q[9], q[10];
ccx q[0], q[8], q[9];
ccx q[1], q[9], q[10];
ccx q[2], q[10], q[11];
cx q[11], q[3];
ccx q[2], q[10], q[11];
cx q[10], q[2];
ccx q[1], q[9], q[10];
cx q[9], q[1];
ccx q[0], q[8], q[9];
cx q[8], q[0];
cx q[9], q[10];
cx q[10], q[11];
cx q[9], q[1];
cx q[10], q[2];
cx q[11], q[3];
x q[0];
x q[1];
x q[2];
x q[3];
cx q[12], q[8];
cx q[13], q[9];
cx q[14], q[10];
cx q[15], q[11];