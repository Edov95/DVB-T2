function [u_hat, iteration] = decode(r,H, N_ldpc, rate, sigmaw2)
%DECODE Apply de message passing algorithm
%   Detailed explanation goes here
    K_ldpc = size(H,2) - size(H,1);
% Initialize the alogrith
    % LLR is still the same since the probability of getting 0 and 1 are
    % the same
    nIterationMax = 50;
    L = zeros(1,N_ldpc); % Each element is the LLR of the corresponding variable node
    M = ones(N_ldpc-K_ldpc,1) * (-2*r/sigmaw2)';   % Initialize the LLR
    LLR = zeros(N_ldpc - K_ldpc, N_ldpc);
    VN = cell(1,N_ldpc);    % Element i contains the indexes of the check nodes in which the variable i is involved
    CN = cell(1,N_ldpc - K_ldpc);  % Element i contains the indexes of the variable nodes involved in the check i
    for vn=1:N_ldpc
        VN{vn} = find(H(:,vn))';
    end
    for cn=1:(N_ldpc - K_ldpc)
        CN{cn} = find(H(cn,:));
    end 
    
% Decoding algo
    for i=1:nIterationMax
        iteration = i;
        disp(i);
        % For each check node compute the messages to the variable nodes        
        for cn=1:(N_ldpc-K_ldpc)
            for vn=CN{cn}                   
                inM = M(cn,CN{cn}(CN{cn}~=vn)); % Messages from all the relevant variables nodes except vn  
                temp = lntanh(abs(inM));
                LLR(cn,vn) = prod(sign(inM))*(lntanh(sum(temp)));                
            end            
        end 
        % For each variable nodes compute the LLR as a sum of intrinsinc and extrinsing information
        for vn=1:N_ldpc
            L(vn) = sum(LLR(VN{vn},vn)) - 2*r(vn)/sigmaw2;             
        end
        
        % Get the estimated output from the llr (Hard decision)
        yCap = zeros(N_ldpc,1);
        yCap(L<0) = 1;
        
        % If the estimated codeword satisfies
        % all the check rules break the cycle        
        if(sum(mod(H*yCap,2)) == 0)
            disp("ok");
            break;
        else                                   
            for vn=1:N_ldpc
                for cn=VN{vn}                    
                    M(cn,vn) = sum(LLR(VN{vn}(VN{vn}~=cn),vn)) - 2*r(vn)/sigmaw2;
                end
            end 
        end
    end
    
    u_hat = yCap;

end

