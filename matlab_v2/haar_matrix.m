function [H] = haar_matrix(n)
%HAAR_MATRIX Summary of this function goes here
%   Detailed explanation goes here

% check input parameter and make sure it's the power of 2
levels = fix(log2(n));
if 2^levels ~= n
    error('please ensure the value of input parameter is a power of 2');

end 
%Initialization
H = 1;
NC = 1/sqrt(2);%normalization constant
LP = [1 1]; 
HP = [1 -1];
% iteration from H=[1] 
for i=1:levels
    H=NC*[kron(H,LP);kron(eye(size(H)),HP)];
end

end

