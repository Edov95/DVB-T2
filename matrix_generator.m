% Code that construct the generating matrix of the LDPC code

clear all; close all; clc;

rate = 37/45;
N_ldpc = 16200;
K_ldpc = N_ldpc*rate;

x = rand(K_ldpc,1);
x(find(x >= 0.5)) = 1;
x(find(x < 0.5)) = 0;

B = generate_B_matrix(rate,N_ldpc);

parity = B*x;
parity = cumsum(parity);
parity = mod(parity,2);

u = [x' parity']';

save('B_56.mat','B');

