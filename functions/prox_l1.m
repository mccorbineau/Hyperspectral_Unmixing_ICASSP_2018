function p = prox_l1(X,gamma) 

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  x: point at which we want to compute the proximity operator
%         gamma: coefficient
%
% Output: p: proximity operator of the l1 norm at X
%
% This function computes the proximity operator of gamma*l1 norm (soft thresholding)
%====================================================================

p = (abs(X)>gamma).*(X-gamma*sign(X));

end