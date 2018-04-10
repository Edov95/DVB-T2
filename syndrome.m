function [csi] = syndrome(r,B, N_ldpc)
%DECODE Calculate the syndrome of the encoded word

% Calc the C Matrix for generate the A matrix
% [null K_ldpc] = size(B);
% C = circshift(speye(N_ldpc - K_ldpc),1) + speye(N_ldpc - K_ldpc);
% C(1,N_ldpc - K_ldpc) = 0;
% % C = abs(mod(inv(C),2));
% % 
% % A = C * B;
% % A = mod(A,2);
% 
% H = [B C];

H = t2_std_dvbt2_bl_ldpc(37/45,N_ldpc,'1.0.1');

csi = H * r;

csi = mod(csi,2);

end

