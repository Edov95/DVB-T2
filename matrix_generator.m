% Code that construct the generating matrix of the LDPC code

clear all; close all;

% Open the related file
fid = fopen('./codes/64800/rate_12');

% Create the systematic matrix
G = sparse(64800,34200);

% Add the identity in the first block
G = sparse(1:34200,1:34200,1);

% Block counter
n_block = 0;

% Populate the matrix
while ~feof(fid)
     x = sscanf(fgetl(fid),'%f');
     for i = 1:360
         G(34200 + x + 1, n_block*360 + i) = 1;
         x = mod(x + i*90, 34200);
         %disp(i.');
     end
     n_block = n_block + 1;
end

fclose(fid);

