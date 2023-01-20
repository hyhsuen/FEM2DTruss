% Created with:       FlExtract v1.13
% Element type:       truss
% Number of nodes:    19
% Number of elements: 34

clear all;

% Node coordinates: x, y
X = [
0	0
0	0.5
0.5	0
0.5	0.5
0.5	1
1	0
1	0.5
1.5	0
1.5	0.5
2	0
2	0.5
2.5	0
2.5	0.5
3	0
3	0.5
3.5	0
3.5	0.5
4	0
4	0.5
];
% Element connectivity: node1_id, node2_id, material_id
IX = [
2	1	1
4	1	1
3	1	1
4	2	1
4	3	1
7	3	1
6	3	1
5	4	1
7	4	1
7	6	1
9	6	1
8	6	1
9	7	1
9	8	1
11	8	1
10	8	1
11	9	1
11	10	1
13	10	1
12	10	1
13	11	1
13	12	1
15	12	1
14	12	1
15	13	1
15	14	1
17	14	1
16	14	1
17	15	1
17	16	1
19	16	1
18	16	1
19	17	1
19	18	1
];
% Element properties: Young's modulus, area
mprop = [
7e+010	0.0002
];
% Nodal diplacements: node_id, degree of freedom (1 - x, 2 - y), displacement
bound = [
2	1	0
2	2	0
5	1	0
5	2	0
];
% Nodal loads: node_id, degree of freedom (1 - x, 2 - y), load
loads = [
18	2	-10000
];
% Control parameters
plotdof = 2;
