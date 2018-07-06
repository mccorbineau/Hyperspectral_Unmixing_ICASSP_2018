function WX = Wav_mult(Wav,X)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  X: variable with multiple channels/endmembers
%         Wav: wavelet transform operator
%
% Output: WX: wavelet coefficients of X
%
% This function computes the wavelet coefficients for every channel of X.
%====================================================================

P  = size(X,1);      % number of endmembers
WX = zeros(size(X)); 

for p = 1:P 
    WX(p,:) = Wav(X(p,:)); % wavelet transform for each channel
end

end
