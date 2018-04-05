clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% 0 SET THE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%

% The steps 1, 2, 3 and 4 must be inserted on a for cycle


% CODES PART

R = 2/3; % The codes rates
% load_matricies(R);

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

% SIGNAL PART

EbN0dBPass = 0.25;          % The pass for the bit SNR
EbN0dB = 0.5:EbN0dBPass:2;  % The bit SNR test values

SNRdB = EbN0dB + 10*log10(2*R); %The SNR values for the signal

LLR = zeros(1,N_ldpc);

% foreach code
%   foreach SNR

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
r = awgn(s,SNRdB(1));

% diff = r - s;
% size(find(diff ~= 0))

%%%%%%%%%%%%
% 4 DECODE %
%%%%%%%%%%%%

% Sono valori decrescenti ed è giusto così perché l'SNR aumenta ma la
% potenza del segnale rimane costante quindi può solo che diminuire la
% potenza dell'errore
sigma2 = 10.^(1./(10*SNRdB(1)));

% Inizialize the LLR for soft decoding
LLR = -2*r/sigma2;

% decode the codeword
u_hat = decoder(r, LLR, B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 COMPARE THE RESULTS AND PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



