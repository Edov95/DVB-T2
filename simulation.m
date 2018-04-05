clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% 0 SET THE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%

% The steps 1, 2, 3 and 4 must be inserted on a for cycle


% CODES PART

R = 37/45; % The codes rates
% load_matricies(R);

load('B_56.mat','B');
[rows,cols] = size(B);

%length of the coded word
N_ldpc = 16200;

% K_ldpc = cols
% Since K rapresent the length of the uncoded word, the numbers of the
% columns of the matrix and the K are the same
K_ldpc = cols;

% Inizialization of the encoded and recived words
d = zeros(N_ldpc,1);
s = zeros(N_ldpc,1);
r = zeros(N_ldpc,1);

% Inizialization of the word
u = zeros(K_ldpc,1);

parity_bits = zeros(N_ldpc - K_ldpc,1);

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
% parity_bits = C * parity_bits;
parity_bits = mod(parity_bits,2);

d(1:K_ldpc) = u;
d(K_ldpc + 1:N_ldpc) = parity_bits;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 SIMULATE THE CHANNEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Map the bits into the signal constellation
% s = bpsk(d);

% Add the noise
% r = awgn(s,SNRdB(1));
r = bsc(s,0.005);

diff = r - s;
size(find(diff ~= 0))

%%%%%%%%%%%%
% 4 DECODE %
%%%%%%%%%%%%

% calculate the syndrome
csi = syndrome(r, B, N_ldpc);

% Sono valori decrescenti ed è giusto così perché l'SNR aumenta ma la
% potenza del segnale rimane costante quindi può solo che diminuire la
% potenza dell'errore
sigma2 = 10.^(1./(10*SNRdB(1)));

% Inizialize the LLR for soft decoding
LLR = -2*r/sigma2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 COMPARE THE RESULTS AND PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



