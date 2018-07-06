function X = Wav_inv_mult(Wav_inv,WX)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  WX: wavelet coefficients of variable X with multiple channels
%         Wav_inv: inverse wavelet transform operator
%
% Output: X: variable whose wavelet coefss are WX
%
% This function applies the inverse wavelet transforms to every channel of WX.
%====================================================================

P = size(WX,1);      % number of channels
X = zeros(size(WX));

for p = 1:P 
    X(p,:) = Wav_inv(WX(p,:)); % inverse wavelet transform for each channel
end
