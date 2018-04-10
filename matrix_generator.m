% Code that construct the generating matrix of the LDPC code

clear all; close all; clc;

rate = 11/15;
R = rate;
N_ldpc = 16200;
K_ldpc = N_ldpc*rate;

x = rand(K_ldpc,1);
x(find(x >= 0.5)) = 1;
x(find(x < 0.5)) = 0;

B = generate_B_matrix(rate,N_ldpc);

save('B_34.mat','B','R');

