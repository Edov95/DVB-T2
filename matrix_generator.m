% Code that construct the generating matrix of the LDPC code

clear all; close all; clc;

rate = 2/3;
N_ldpc = 64800;
K_ldpc = 64800*rate;

x = rand(K_ldpc,1);
x(find(x >= 0.5)) = 1;
x(find(x < 0.5)) = 0;

B = generate_B_matrix(rate,N_ldpc);

parity = B*x;
parity = cumsum(parity);
parity = mod(parity,2);



