function [csi] = syndrome(r,B, N_ldpc)
%DECODE Calculate the syndrome of the encoded word

% Calc the C Matrix for generate the A matrix
[null K_ldpc] = size(B);
C = circshift(speye(N_ldpc - K_ldpc),1) + speye(N_ldpc - K_ldpc);
C(1,N_ldpc - K_ldpc) = 0;

A = C * B;

H = [A speye(N_ldpc-K_ldpc)];

csi = H * r;

csi = mod(csi,2);

end

