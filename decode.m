function [u_hat] = decode(r,LLR,B, N_ldpc, rate)
%DECODE Apply de message passing algorithm
%   Detailed explanation goes here
K_ldpc = N_ldpc * rate;
C = circshift(speye(N_ldpc - K_ldpc),1) + speye(N_ldpc - K_ldpc);
C(1,N_ldpc - K_ldpc) = 0;

A = C * B;

H = [A speye(N_ldpc-K_ldpc)];



end

