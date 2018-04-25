clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%
% 0 SET THE PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%

%length of the coded word
N_ldpc = 16200;

% Inizialization of the encoded and recived words
d = zeros(N_ldpc,1);
s = zeros(N_ldpc,1);
r = zeros(N_ldpc,1);


Nit = 1e5; % iterations number: try with Nit diff numbers of codewords

%Riscrivere meglio (forse inutile) e giocare con le matrici sparse
%     [max_check_degree,check_node_ones,max_variable_degree,variable_node_ones]=one_finder(H);

% Inizialization of the word
u_hat = zeros(N_ldpc,1);

% SIGNAL PART
EbN0dBPass = 1;          % The pass for the bit SNR
EbN0dB = 6:EbN0dBPass:14;  % The bit SNR test values

LLR = zeros(1,N_ldpc);

npack = zeros([length(EbN0dB),4]);
%     nerr_2 = zeros(size(EbN0dB));
%     npack_2 = zeros(size(EbN0dB));

%%%%%%%%%%%%%%
% SIMULATION %
%%%%%%%%%%%%%%

tic

for r_i = 1:4
    
    switch r_i
        case 1 
            load('B_12.mat','B','R');
        case 2
            load('B_23.mat','B','R');
        case 3
            load('B_34.mat','B','R');
        case 4
            load('B_56.mat','B','R');
    end
    
    K_ldpc = R * N_ldpc;
    
    prec = 1+log10(1/(K_ldpc*Nit)); % expected precision
    disp(['Expected precision = 10^' num2str(prec)])

    C = circshift(speye(N_ldpc - K_ldpc),1) + speye(N_ldpc - K_ldpc);
    C(1,N_ldpc - K_ldpc) = 0;

    H = [B C];
    
    SNRdB = EbN0dB + 10*log10(2*R); %The SNR values for the signal
    
    u = zeros(K_ldpc,1);
    
    ldpcFEC = comm.LDPCDecoder('ParityCheckMatrix',H, 'MaximumIterationCount', 50, 'IterationTerminationCondition', 'Parity check satisfied', ...
    'DecisionMethod', 'Hard decision', 'NumIterationsOutputPort', true);

    parity_bits = zeros(N_ldpc - K_ldpc,1);
    
    for it = 1:Nit
        fprintf("Iteration :%i\n",it);
        
        % 1 GENERATE A CODE WORD
        
        % The uncoded word
        u = randi(2,K_ldpc,1) - 1;
        
        % 2 ENCODE
        parity_bits = mod(cumsum(B * u),2);
        
        d = [u; parity_bits];
        
        % 3 SIMULATE THE CHANNEL
        
        % Map the bits into the signal constellation
        s = bpsk(d);
        
        for m = 1:length(SNRdB)
            
            % Add the noise
            r = awgn(s,SNRdB(m));
            
            sigma2 = 10.^(-SNRdB(m)/10);
            
            % Inizialize the LLR for soft decoding
            LLR = -2*r/sigma2;
            
            %         [u_hat, iteration] = decode_ldpc_new(200,LLR,check_node_ones,max_check_degree,variable_node_ones,max_variable_degree,N_ldpc - K_ldpc,N_ldpc);
            %         error_u_hat = sum(u_hat(1:K_ldpc)~=u);
            %         if(error_u_hat > 0)
            %             npack(m) = npack(m) + 1;
            %         end
            %
            [u_hat, iteration2] = step(ldpcFEC, LLR);
            if(sum(u_hat ~= u) > 0)
                npack(m,r_i) = npack(m,r_i) + 1;
            end
            %         fprintf("Error in u_hat_2 %i\n",error_u_hat_2);
            %         fprintf("Corrected %i errors with MATLAB decoder\n",error_r_hat - error_u_hat_2);
        end
    end
    clear B R C H;
end
toc

fer = npack./Nit;
LDPCplot(EbN0dB,fer);


