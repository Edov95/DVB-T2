function [ y ] = lntanh( x )
%LNTANH 
 
 y = -log(tanh(x/2));
 y(y==Inf) = 11;
 y(find(y>10))=3;
 y(find(x>10))=0;
 
end

