function [B] = generate_B_matrix(rate,N_ldpc)
%GENERATE_B_MATRIX Summary of this function goes here
%   Given the code rate and the lenght of the word return the matrix B
%   for that code

    K_ldpc = N_ldpc * rate; % length of the information word
    parity_length = N_ldpc - K_ldpc; % length of the parity bits


    % Open the related file
    % need to use a switch condition to choose the correct code
    fid = fopen('./codes/rate_34');

    % Create the systematic matrix
    B = sparse(N_ldpc-K_ldpc,K_ldpc);

    % Block counter
    n_block = 0;

    % Populate the matrix
    while ~feof(fid)
         x = sscanf(fgetl(fid),'%f');
         for i = 1:360
             B(x + 1, n_block*360 + i) = 1;
             x = mod(x + i*12, parity_length);
             %disp(i.');
         end
         n_block = n_block + 1;
    end

    fclose(fid);

end

