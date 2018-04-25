% Code that construct the generating matrix of the LDPC code

clear all; close all; clc;

rate = 2/3;
R = rate;
N_ldpc = 16200;

B = generate_B_matrix(rate,N_ldpc);

save('B_23.mat','B','R');

