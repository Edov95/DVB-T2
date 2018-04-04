clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% 0 SET THE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%

% The steps 1, 2, 3 and 4 must be inserted on a for cycle

load('B_23.mat','B');
[rows,cols] = size(B);

%length of the coded word
N_ldpc = 64800;

% K_ldpc = cols
% Since K rapresent the length of the uncoded word, the numbers of the
% columns of the matrix and the K are the same
K_ldpc = cols;

% Inizialization of the encoded and recived words
d = zeros(1,N_ldpc);
s = zeros(1,N_ldpc);
r = zeros(1,N_ldpc);

% Inizialization of the word
u = zeros(K_ldpc,1);

parity_bits = zeros(1,N_ldpc - K_ldpc);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 GENERATE A CODE WORD %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% The uncoded word
u = rand(K_ldpc,1);
u(find(u >= 0.5)) = 1;
u(find(u < 0.5)) = 0;

%%%%%%%%%%%%
% 2 ENCODE %
%%%%%%%%%%%%

parity_bits = B * u;
parity_bits = cumsum(parity_bits);
parity_bits = mod(parity_bits,2);

d(1:K_ldpc) = u;
d(K_ldpc + 1:N_ldpc) = parity_bits;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 SIMULATE THE CHANNEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Map the bits into the signal constellation
s = bpsk(d);

% Add the channel noise
r = send_over_channel(s);

%%%%%%%%%%%%
% 4 DECODE %
%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 COMPARE THE RESULTS AND PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%