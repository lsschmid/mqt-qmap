// Benchmark was created by MQT Bench on 2023-06-29
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: v1.0.0
// Qiskit version: {'qiskit-terra': '0.24.1', 'qiskit-aer': '0.12.0', 'qiskit-ignis': None, 'qiskit-ibmq-provider': '0.20.2', 'qiskit': '0.43.1', 'qiskit-nature': '0.6.2', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.5.0', 'qiskit-machine-learning': '0.6.1'}
// Used Gate Set: ['rx', 'rz', 'cz', 'measure']

OPENQASM 2.0;
include "qelib1.inc";
qreg q[50];
creg meas[50];
rz(-pi/2) q[0];
rx(-pi/4) q[0];
rz(-pi/2) q[1];
rx(-0.955316618124509) q[1];
rz(-pi/2) q[2];
rx(-pi/3) q[2];
rz(-pi/2) q[3];
rx(-1.1071487177940904) q[3];
rz(-pi/2) q[4];
rx(-1.1502619915109318) q[4];
rz(-pi/2) q[5];
rx(-1.1831996401397162) q[5];
rz(-pi/2) q[6];
rx(-1.209429202888189) q[6];
rz(-pi/2) q[7];
rx(-1.2309594173407747) q[7];
rz(-pi/2) q[8];
rx(-1.2490457723982544) q[8];
rz(-pi/2) q[9];
rx(-1.2645189576252271) q[9];
rz(-pi/2) q[10];
rx(-1.277953555066321) q[10];
rz(-pi/2) q[11];
rx(-1.2897614252920828) q[11];
rz(-pi/2) q[12];
rx(-1.3002465638163239) q[12];
rz(-pi/2) q[13];
rx(-1.3096389158918722) q[13];
rz(-pi/2) q[14];
rx(-1.3181160716528182) q[14];
rz(-pi/2) q[15];
rx(-1.3258176636680323) q[15];
rz(-pi/2) q[16];
rx(-1.3328552019646882) q[16];
rz(-pi/2) q[17];
rx(-1.3393189628247182) q[17];
rz(-pi/2) q[18];
rx(-1.3452829208967654) q[18];
rz(-pi/2) q[19];
rx(-1.3508083493994372) q[19];
rz(-pi/2) q[20];
rx(-1.3559464937191843) q[20];
rz(-pi/2) q[21];
rx(-1.3607405877236578) q[21];
rz(-pi/2) q[22];
rx(-1.3652273956337226) q[22];
rz(-pi/2) q[23];
rx(-1.369438406004566) q[23];
rz(-pi/2) q[24];
rx(-1.3734007669450157) q[24];
rz(-pi/2) q[25];
rx(-1.37713802635057) q[25];
rz(-pi/2) q[26];
rx(-1.3806707234484301) q[26];
rz(-pi/2) q[27];
rx(-1.3840168657133032) q[27];
rz(-pi/2) q[28];
rx(-1.3871923165159783) q[28];
rz(-pi/2) q[29];
rx(-1.3902111126041987) q[29];
rz(-pi/2) q[30];
rx(-1.393085725949785) q[30];
rz(-pi/2) q[31];
rx(-1.3958272811292074) q[31];
rz(-pi/2) q[32];
rx(-1.3984457368955736) q[32];
rz(-pi/2) q[33];
rx(-1.400950038711223) q[33];
rz(-pi/2) q[34];
rx(-1.4033482475752073) q[34];
rz(-pi/2) q[35];
rx(-1.4056476493802696) q[35];
rz(-pi/2) q[36];
rx(-1.4078548481843771) q[36];
rz(-pi/2) q[37];
rx(-1.4099758461204321) q[37];
rz(-pi/2) q[38];
rx(-1.4120161121491361) q[38];
rz(-pi/2) q[39];
rx(-1.4139806414504958) q[39];
rz(-pi/2) q[40];
rx(-1.4158740069240832) q[40];
rz(-pi/2) q[41];
rx(-1.4177004040080423) q[41];
rz(-pi/2) q[42];
rx(-1.419463689817681) q[42];
rz(-pi/2) q[43];
rx(-1.4211674174353792) q[43];
rz(-pi/2) q[44];
rx(-1.422814866046113) q[44];
rz(-pi/2) q[45];
rx(-1.4244090675006476) q[45];
rz(-pi/2) q[46];
rx(-1.425952829796337) q[46];
rz(-pi/2) q[47];
rx(-1.4274487578895314) q[47];
rz(-pi/2) q[48];
rx(-1.4288992721907328) q[48];
rx(pi) q[49];
cz q[49],q[48];
rx(1.4288992721907323) q[48];
cz q[48],q[47];
rx(1.4274487578895305) q[47];
cz q[47],q[46];
rx(1.4259528297963369) q[46];
cz q[46],q[45];
rx(1.4244090675006473) q[45];
cz q[45],q[44];
rx(1.4228148660461128) q[44];
cz q[44],q[43];
rx(1.421167417435379) q[43];
cz q[43],q[42];
rx(1.4194636898176807) q[42];
cz q[42],q[41];
rx(1.417700404008042) q[41];
cz q[41],q[40];
rx(1.415874006924083) q[40];
cz q[40],q[39];
rx(1.4139806414504958) q[39];
cz q[39],q[38];
rx(1.4120161121491355) q[38];
cz q[38],q[37];
rx(1.4099758461204315) q[37];
cz q[37],q[36];
rx(1.4078548481843771) q[36];
cz q[36],q[35];
rx(1.4056476493802692) q[35];
cz q[35],q[34];
rx(1.4033482475752075) q[34];
cz q[34],q[33];
rx(1.4009500387112224) q[33];
cz q[33],q[32];
rx(1.3984457368955732) q[32];
cz q[32],q[31];
rx(1.3958272811292074) q[31];
cz q[31],q[30];
rx(1.3930857259497846) q[30];
cz q[30],q[29];
rx(1.3902111126041987) q[29];
cz q[29],q[28];
rx(1.3871923165159783) q[28];
cz q[28],q[27];
rx(1.3840168657133025) q[27];
cz q[27],q[26];
rx(1.3806707234484301) q[26];
cz q[26],q[25];
rx(1.37713802635057) q[25];
cz q[25],q[24];
rx(1.373400766945016) q[24];
cz q[24],q[23];
rx(1.3694384060045655) q[23];
cz q[23],q[22];
rx(1.3652273956337222) q[22];
cz q[22],q[21];
rx(1.3607405877236578) q[21];
cz q[21],q[20];
rx(1.3559464937191839) q[20];
cz q[20],q[19];
rx(1.3508083493994372) q[19];
cz q[19],q[18];
rx(1.345282920896765) q[18];
cz q[18],q[17];
rx(1.339318962824718) q[17];
cz q[17],q[16];
rx(1.3328552019646878) q[16];
cz q[16],q[15];
rx(1.3258176636680323) q[15];
cz q[15],q[14];
rx(1.3181160716528175) q[14];
cz q[14],q[13];
rx(1.309638915891872) q[13];
cz q[13],q[12];
rx(1.3002465638163239) q[12];
cz q[12],q[11];
rx(1.2897614252920826) q[11];
cz q[11],q[10];
rx(1.2779535550663208) q[10];
cz q[10],q[9];
rz(pi/2) q[49];
rx(pi/2) q[49];
cz q[48],q[49];
rx(-pi/2) q[48];
rz(-pi) q[48];
cz q[47],q[48];
rx(-pi/2) q[47];
rz(-pi) q[47];
cz q[46],q[47];
rx(-pi/2) q[46];
rz(-pi) q[46];
cz q[45],q[46];
rx(-pi/2) q[45];
rz(-pi) q[45];
cz q[44],q[45];
rx(-pi/2) q[44];
rz(-pi) q[44];
cz q[43],q[44];
rx(-pi/2) q[43];
rz(-pi) q[43];
cz q[42],q[43];
rx(-pi/2) q[42];
rz(-pi) q[42];
cz q[41],q[42];
rx(-pi/2) q[41];
rz(-pi) q[41];
cz q[40],q[41];
rx(-pi/2) q[40];
rz(-pi) q[40];
cz q[39],q[40];
rx(-pi/2) q[39];
rz(-pi) q[39];
cz q[38],q[39];
rx(-pi/2) q[38];
rz(-pi) q[38];
cz q[37],q[38];
rx(-pi/2) q[37];
rz(-pi) q[37];
cz q[36],q[37];
rx(-pi/2) q[36];
rz(-pi) q[36];
cz q[35],q[36];
rx(-pi/2) q[35];
rz(-pi) q[35];
cz q[34],q[35];
rx(-pi/2) q[34];
rz(-pi) q[34];
cz q[33],q[34];
rx(-pi/2) q[33];
rz(-pi) q[33];
cz q[32],q[33];
rx(-pi/2) q[32];
rz(-pi) q[32];
cz q[31],q[32];
rx(-pi/2) q[31];
rz(-pi) q[31];
cz q[30],q[31];
rx(-pi/2) q[30];
rz(-pi) q[30];
cz q[29],q[30];
rx(-pi/2) q[29];
rz(-pi) q[29];
cz q[28],q[29];
rx(-pi/2) q[28];
rz(-pi) q[28];
cz q[27],q[28];
rx(-pi/2) q[27];
rz(-pi) q[27];
cz q[26],q[27];
rx(-pi/2) q[26];
rz(-pi) q[26];
cz q[25],q[26];
rx(-pi/2) q[25];
rz(-pi) q[25];
cz q[24],q[25];
rx(-pi/2) q[24];
rz(-pi) q[24];
cz q[23],q[24];
rx(-pi/2) q[23];
rz(-pi) q[23];
cz q[22],q[23];
rx(-pi/2) q[22];
rz(-pi) q[22];
cz q[21],q[22];
rx(-pi/2) q[21];
rz(-pi) q[21];
cz q[20],q[21];
rx(-pi/2) q[20];
rz(-pi) q[20];
cz q[19],q[20];
rx(-pi/2) q[19];
rz(-pi) q[19];
cz q[18],q[19];
rx(-pi/2) q[18];
rz(-pi) q[18];
cz q[17],q[18];
rx(-pi/2) q[17];
rz(-pi) q[17];
cz q[16],q[17];
rx(-pi/2) q[16];
rz(-pi) q[16];
cz q[15],q[16];
rx(-pi/2) q[15];
rz(-pi) q[15];
cz q[14],q[15];
rx(-pi/2) q[14];
rz(-pi) q[14];
cz q[13],q[14];
rx(-pi/2) q[13];
rz(-pi) q[13];
cz q[12],q[13];
rx(-pi/2) q[12];
rz(-pi) q[12];
cz q[11],q[12];
rx(-pi/2) q[11];
rz(-pi) q[11];
cz q[10],q[11];
rx(-pi/2) q[10];
rz(-pi) q[10];
rx(-pi/2) q[11];
rz(-pi/2) q[11];
rx(-pi/2) q[12];
rz(-pi/2) q[12];
rx(-pi/2) q[13];
rz(-pi/2) q[13];
rx(-pi/2) q[14];
rz(-pi/2) q[14];
rx(-pi/2) q[15];
rz(-pi/2) q[15];
rx(-pi/2) q[16];
rz(-pi/2) q[16];
rx(-pi/2) q[17];
rz(-pi/2) q[17];
rx(-pi/2) q[18];
rz(-pi/2) q[18];
rx(-pi/2) q[19];
rz(-pi/2) q[19];
rx(-pi/2) q[20];
rz(-pi/2) q[20];
rx(-pi/2) q[21];
rz(-pi/2) q[21];
rx(-pi/2) q[22];
rz(-pi/2) q[22];
rx(-pi/2) q[23];
rz(-pi/2) q[23];
rx(-pi/2) q[24];
rz(-pi/2) q[24];
rx(-pi/2) q[25];
rz(-pi/2) q[25];
rx(-pi/2) q[26];
rz(-pi/2) q[26];
rx(-pi/2) q[27];
rz(-pi/2) q[27];
rx(-pi/2) q[28];
rz(-pi/2) q[28];
rx(-pi/2) q[29];
rz(-pi/2) q[29];
rx(-pi/2) q[30];
rz(-pi/2) q[30];
rx(-pi/2) q[31];
rz(-pi/2) q[31];
rx(-pi/2) q[32];
rz(-pi/2) q[32];
rx(-pi/2) q[33];
rz(-pi/2) q[33];
rx(-pi/2) q[34];
rz(-pi/2) q[34];
rx(-pi/2) q[35];
rz(-pi/2) q[35];
rx(-pi/2) q[36];
rz(-pi/2) q[36];
rx(-pi/2) q[37];
rz(-pi/2) q[37];
rx(-pi/2) q[38];
rz(-pi/2) q[38];
rx(-pi/2) q[39];
rz(-pi/2) q[39];
rx(-pi/2) q[40];
rz(-pi/2) q[40];
rx(-pi/2) q[41];
rz(-pi/2) q[41];
rx(-pi/2) q[42];
rz(-pi/2) q[42];
rx(-pi/2) q[43];
rz(-pi/2) q[43];
rx(-pi/2) q[44];
rz(-pi/2) q[44];
rx(-pi/2) q[45];
rz(-pi/2) q[45];
rx(-pi/2) q[46];
rz(-pi/2) q[46];
rx(-pi/2) q[47];
rz(-pi/2) q[47];
rx(-pi/2) q[48];
rz(-pi/2) q[48];
rx(-pi/2) q[49];
rz(-pi/2) q[49];
rx(1.264518957625227) q[9];
cz q[9],q[8];
rx(1.2490457723982542) q[8];
cz q[8],q[7];
rx(1.2309594173407745) q[7];
cz q[7],q[6];
rx(1.209429202888189) q[6];
cz q[6],q[5];
rx(1.1831996401397153) q[5];
cz q[5],q[4];
rx(1.1502619915109316) q[4];
cz q[4],q[3];
rx(1.1071487177940902) q[3];
cz q[3],q[2];
rx(pi/3) q[2];
cz q[2],q[1];
rx(0.9553166181245094) q[1];
cz q[1],q[0];
rx(pi/4) q[0];
rz(pi/2) q[0];
cz q[9],q[10];
rx(-pi/2) q[10];
rz(-pi/2) q[10];
rx(-pi/2) q[9];
rz(-pi) q[9];
cz q[8],q[9];
rx(-pi/2) q[8];
rz(-pi) q[8];
cz q[7],q[8];
rx(-pi/2) q[7];
rz(-pi) q[7];
cz q[6],q[7];
rx(-pi/2) q[6];
rz(-pi) q[6];
cz q[5],q[6];
rx(-pi/2) q[5];
rz(-pi) q[5];
cz q[4],q[5];
rx(-pi/2) q[4];
rz(-pi) q[4];
cz q[3],q[4];
rx(-pi/2) q[3];
rz(-pi) q[3];
cz q[2],q[3];
rx(-pi/2) q[2];
rz(-pi) q[2];
cz q[1],q[2];
rx(-pi/2) q[1];
rz(-pi) q[1];
cz q[0],q[1];
rx(-pi/2) q[1];
rz(-pi/2) q[1];
rx(-pi/2) q[2];
rz(-pi/2) q[2];
rx(-pi/2) q[3];
rz(-pi/2) q[3];
rx(-pi/2) q[4];
rz(-pi/2) q[4];
rx(-pi/2) q[5];
rz(-pi/2) q[5];
rx(-pi/2) q[6];
rz(-pi/2) q[6];
rx(-pi/2) q[7];
rz(-pi/2) q[7];
rx(-pi/2) q[8];
rz(-pi/2) q[8];
rx(-pi/2) q[9];
rz(-pi/2) q[9];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15],q[16],q[17],q[18],q[19],q[20],q[21],q[22],q[23],q[24],q[25],q[26],q[27],q[28],q[29],q[30],q[31],q[32],q[33],q[34],q[35],q[36],q[37],q[38],q[39],q[40],q[41],q[42],q[43],q[44],q[45],q[46],q[47],q[48],q[49];
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
measure q[10] -> meas[10];
measure q[11] -> meas[11];
measure q[12] -> meas[12];
measure q[13] -> meas[13];
measure q[14] -> meas[14];
measure q[15] -> meas[15];
measure q[16] -> meas[16];
measure q[17] -> meas[17];
measure q[18] -> meas[18];
measure q[19] -> meas[19];
measure q[20] -> meas[20];
measure q[21] -> meas[21];
measure q[22] -> meas[22];
measure q[23] -> meas[23];
measure q[24] -> meas[24];
measure q[25] -> meas[25];
measure q[26] -> meas[26];
measure q[27] -> meas[27];
measure q[28] -> meas[28];
measure q[29] -> meas[29];
measure q[30] -> meas[30];
measure q[31] -> meas[31];
measure q[32] -> meas[32];
measure q[33] -> meas[33];
measure q[34] -> meas[34];
measure q[35] -> meas[35];
measure q[36] -> meas[36];
measure q[37] -> meas[37];
measure q[38] -> meas[38];
measure q[39] -> meas[39];
measure q[40] -> meas[40];
measure q[41] -> meas[41];
measure q[42] -> meas[42];
measure q[43] -> meas[43];
measure q[44] -> meas[44];
measure q[45] -> meas[45];
measure q[46] -> meas[46];
measure q[47] -> meas[47];
measure q[48] -> meas[48];
measure q[49] -> meas[49];