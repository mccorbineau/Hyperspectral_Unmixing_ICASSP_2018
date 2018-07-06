function [approx, detail] = sort_wavelet_coeffs(Wx,nRow,nCol,P,J)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  Wx: wavelet coefficients
%         nRow, nCol: number of rows and columns, resp.
%         P: number of endmembers (or channels)
%         J: decomposition level
%
% Output: approx: approximation wavelet coefficients
%         detail: detail wavelet coefficients
%
% This function separates the detail coeffs from the approx coeffs of 
% the wavelet decomposition.
%====================================================================

Wx           = reshape(Wx',nRow,nCol,P);
nRow_approx  = nRow/(2^J);
nCol_approx  = nCol/(2^J);

approx                                 =  Wx(1:nRow_approx,1:nCol_approx,:);
detail                                 =  Wx;
detail(1:nRow_approx,1:nCol_approx,:)  =  zeros(nRow_approx,nCol_approx,P);

end

