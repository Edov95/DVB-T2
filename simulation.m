clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% 0 SET THE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%

% The steps 1, 2, 3 and 4 must be inserted on a for cycle

% DEBUG
	% set to 1 for debugging and to 0 for linear execution
    DEBUG = 0;

% CODES PART
    R = 37/45; % The codes rates
    
    load('B_56.mat','B');

    %length of the coded word
    N_ldpc = 16200;
    
        % Inizialization of the encoded and recived words
    d = zeros(N_ldpc,1);
    s = zeros(N_ldpc,1);
    r = zeros(N_ldpc,1);
    
    K_ldpc = R * N_ldpc;
    
    Nit = 1e2; % iterations number: try with Nit diff numbers of codewords
    
    prec = 1+log10(1/(K_ldpc*Nit)); % expected precision
    disp(['Expected precision = 10^' num2str(prec)])
    
    C = circshift(speye(N_ldpc - K_ldpc),1) + speye(N_ldpc - K_ldpc);
    C(1,N_ldpc - K_ldpc) = 0;
    
    H = [B C];
    
    %Riscrivere meglio (forse inutile) e giocare con le matrici sparse
%     [max_check_degree,check_node_ones,BIGVALUE_COLS,max_variable_degree,variable_node_ones,BIGVALUE_ROWS]=one_finder(H);
    
%     Posso usarlo?
    ldpcDEC = comm.LDPCDecoder('ParityCheckMatrix',H, 'MaximumIterationCount', 100, 'IterationTerminationCondition', 'Parity check satisfied', ...
    'DecisionMethod', 'Soft decision', 'NumIterationsOutputPort', true, 'FinalParityChecksOutputPort', true);
    
    % Inizialization of the word
    u = zeros(K_ldpc,1);
    u_hat = zeros(N_ldpc,1);

    parity_bits = zeros(N_ldpc - K_ldpc,1);

% SIGNAL PART
    EbN0dBPass = 0.25;          % The pass for the bit SNR
    EbN0dB = 7:EbN0dBPass:11;  % The bit SNR test values

    SNRdB = EbN0dB + 10*log10(2*R); %The SNR values for the signal

    LLR = zeros(1,N_ldpc);
    
    nerr = zeros(size(EbN0dB));
    npack = zeros(size(EbN0dB));

%%%%%%%%%%%%%%   
% SIMULATION %
%%%%%%%%%%%%%%

tic

for it = 1:Nit
    fprintf("Iteration :%i\n",it);
    
    % 1 GENERATE A CODE WORD

    % The uncoded word
    u = rand(K_ldpc,1);
    u(u >= 0.5) = 1;
    u(u < 0.5) = 0;

    % 2 ENCODE    
    parity_bits = mod(cumsum(B * u),2);

    d = [u; parity_bits];

    for m = 1:length(SNRdB)
            
        % 3 SIMULATE THE CHANNEL
        
        % Map the bits into the signal constellation
        s = bpsk(d);
        
        % Add the noise
        r = awgn(s,SNRdB(m));
        
        % 4 DECODE
        
        % Sono valori decrescenti ed è giusto così perché l'SNR aumenta ma la
        % potenza del segnale rimane costante quindi può solo che diminuire la
        % potenza dell'errore
        sigma2 = 10.^(-(SNRdB(m)/10));
        
        % Inizialize the LLR for soft decoding
        LLR = -2*r/sigma2;
        
        %             [u_hat,iteration]=decode_ldpc_new(200,LLR,check_node_ones,max_check_degree,BIGVALUE_COLS-1,variable_node_ones,max_variable_degree,BIGVALUE_ROWS-1,N_ldpc - K_ldpc,N_ldpc);
        %             u_hat = transpose(u_hat);
        [u_hat, iteration, finalParityChecks] = step(ldpcDEC, LLR);
        u_hat = 1 * (u_hat <=0);
        
        % count errors
        nerr(m) = nerr(m) + sum(u_hat(1:K_ldpc)~=u);
        npack(m) = npack(m) + 1;
      
    end
end
toc

% 5 COMPARE THE RESULTS AND PLOTS
Pbit = nerr./(npack*K_ldpc);
% expected BER
Gamma = 10.^(SNRdB/10); % SNR
Q = @(x) 0.5*erfc((x)/sqrt(2));
Pexp = Q(sqrt(5*Gamma))+4*Q(sqrt(6*Gamma))+12*Q(sqrt(7*Gamma));
% show results
figure;
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
set(gca,'FontSize',14);
semilogy(SNRdB,Pbit,'k-',...
    SNRdB,Pexp,'g--',...
    SNRdB,ones(size(SNRdB))*10^prec,'k-');
axis([min(SNRdB) max(SNRdB) 1e-8 1e0])
legend('Montecarlo','Closed form')
xlabel('SNR $\Gamma$  [dB]')
ylabel('BER $P_{\rm bit}$')
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on',...
    'YGrid', 'on', 'XGrid', 'on');


