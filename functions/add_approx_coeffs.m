function [ Wx ] = add_approx_coeffs( P,J,nRow,nCol,approx,detail )

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  P: number of endmembers (or materials, or channels)
%         J: decomposition level   
%         nRow, nCol: image is of size nRow x nCol
%         approx: approximation coefficients of wavelet decomposition 
%                 (size n_Row/(2^J) x n_Col/(2^J) x P)
%         detail: detail coefficients of wavelet decomposition 
%                 (size nRow x nCol x P)
%
% From separate approximation and detail coefficients, this function 
% re-builds the full wavelet decomposition
%====================================================================

nRow_approx                            =  nRow/(2^J);
nCol_approx                            =  nCol/(2^J);
detail(1:nRow_approx,1:nCol_approx,:)  =  approx;
Wx                                     =  reshape(detail,nRow*nCol,P)';


end

