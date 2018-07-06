function [PX,v] = prox_metric(X,A_1,norm_A_1,gamma,precision,v,nRow,nCol,J,Wav,Wav_inv)

%====================================================================
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  X: point at which the prox must be computed
%         A_1: inverse of the preconditioner   
%         norm_A_1: norm of the inverse of the preconditioning matrix
%         gamma: coeff for the prox
%         precision: precision for computing the prox
%         v: dual variable used for warm-restart
%         nRow, nCol: number of rows and columns, resp.
%         J: decomposition level
%         Wav, Wav_inv: wavelet transform and inverse wavelet transform operators
%
% Output: PX: = argmin_Y 0.5*||X-Y||^2_A + gamma*||Wav*Y||_{1,detail}
%             = argmin_Y 0.5*||A^{1/2}*X-A^{1/2}*Y||^2 + gamma*||Wav*Y||_{1,detail}
%             = A^{-1/2}* argmin_Y' 0.5*||A^{1/2}*X-Y'||^2 + gamma*||Wav*A^{-1/2}*Y'||_{1,detail}
%         v: dual variable used for warm-restart
%
% This function computes the proximity operator of the l1 norm of the
% detail wavelet coefficients, with a preconditioning matrix. 
% The algorithm implemented here is taken from
%
% P.-L. Combettes, D. Dung, and B.C. Vu. Proximity for sums of composite
% functions. Journal of Mathematical Analysis and Applications, vol. 380,
% no. 2, pp. 680-688, 2011. 
%====================================================================

[P,N] = size(X); % number of materials and number of pixels

if(isempty(A_1)) 
    %%%% in case there is no preconditioning matrix
    WX               = Wav_mult(Wav,X);                                % wavelet transform
    [approx,detail]  = sort_wavelet_coeffs(WX,nRow,nCol,P,J);          % identify detail and approx wavelet coeffs
    detail           = prox_l1(detail,gamma);                          % apply the prox of the l1 norm on the detail coeffs
    prox_WX          = add_approx_coeffs(P,J,nRow,nCol,approx,detail); % re-built full wavelet coeffs
    PX               = Wav_inv_mult(Wav_inv,prox_WX);                  % go back to image space
else
    %%%%% start Dual Forward-Backward algorithm
    x = X(:); 
    if (isempty(v))
        
        %%%% if there is no warm-restart
        v = zeros(size(x));
    end
    
    %%%% parameters of DFB
    rho   =  1/(norm_A_1);
    delta =  min(1,rho) - 1e-8;
    step  =  2*rho-delta ;
    
    NbIt = 1000; % maximal number of iterations
    pxold = x;   % store current iterate for the stopping criterion
    
    for i = 1:NbIt
        px = x - A_1(v) ;
        u  = v + step *px;
        v  = u - step .*reshape(prox_metric(reshape(u,P,N)/step,[],[],gamma/step,[],[],nRow,nCol,J,Wav,Wav_inv),P*N,1); 
        if(norm(px - pxold) < precision && i>1); break; end
        pxold = px;
    end
    PX = reshape(px,P,N);
end

 