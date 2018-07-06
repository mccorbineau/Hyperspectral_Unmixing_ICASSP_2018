function [X,obj_vec,snr_vec,X_Xinf_vec,time_vec] = PIPA(Y,S,reg,time_max,nRow,nCol,Xtrue,Xinf)

%====================================================================
% Associated citation:
%
% M.-C. Corbineau, E. Chouzenoux, and J.-C. Pesquet. PIPA: a new proximal
% interior point algorithm for large-scale convex optimization.
% In Proceedings of the 43rd IEEE International Conference on Acoustics, 
% Speech and Signal Processing (ICASSP 2018), Calgary, Canada, 15-20 April 2018.
%
% Kindly report any suggestions or corrections to
% marie-caroline.corbineau@centralesupelec.fr
%
% Input:  Y: data
%         S: library of the spectral signatures of the endmembers
%         reg: regularization parameter
%         time_max: maximal running time
%         nRow, nCol: number of rows and columns, resp.
%         Xtrue: ground-truth of the abundance maps (used only to compute the signal-to-noise ratio)
%         Xinf: solution to the optimization problem (used only to compute the distance to solution)
%
% Output: X: final iterate
%         obj_vec: objective function value for every iteration
%         snr_vec: signal-to-noise ratio at every iteration
%         X_Xinf_vec: normalized distance to solution at every iteration
%         time_vec: time at every iteration
%
% Proximal interior point algorithm for solving:
%  minimize_X      0.5*||S*X-Y||^2_2 + reg*||W*X||_{1,detail}  
%  subject to      for every pixel p, sum_i(X_{i,p})<=1
%                  X >= 0
%
% where W is an orthogonal wavelet transform and ||W*.||_{1,detail} is the 
% l1 norm of the detail coefficients.
%====================================================================

nPix  = nRow*nCol;   % number of pixels
nEnd  = size(S,2);   % number of materials or endmembers

%%% algorithm parameters
maxiter_backtrack  = 1000;  % maximum number of iterations for the backtracking linesearch
gamma_max          = 0.4;   % maximal step size
theta_backtrack    = 0.80;  % granularity of backtracking search
delta_backtrack    = 0.99;  % wolfe parameter of backtracking search
eps1               = 1e5;   % for stopping criterion 1 (for updating the barrier parameter)
eps2               = 1e8;   % for stopping criterion 2
eps3               = 1e8;   % for stopping criterion 3
zeta               = 1;     % used to ensure that stopping criteria decrease faster than barrier parameter
rho                = 1.5;   % geometric decrease for the barrier parameters    
mu                 = 0.01;  % initial barrier parameter
precision          = 1;     % precision for computing the proximity operator
VV                 = [];    % used for warm restart when computing the prox

%%% Wavelet operators
J        = 2;                             % resolution level
qmf      = MakeONFilter('Daubechies',4);  % wavelet 
[~,Jmax] = quadlength(rand(nRow,nCol));
Wav      = @(x) reshape(FWT2_PO(reshape(x,nRow,nCol),Jmax-J,qmf),nPix,1)'; % direct wavelet decomposition operator
Wav_inv  = @(z) reshape(IWT2_PO(reshape(z,nRow,nCol),Jmax-J,qmf),nPix,1)'; % inverse wavelet decomposition operator 

%%% useful variables
C      = [speye(nEnd);-ones(1,nEnd)];     % constraint matrix
c      = [zeros(nEnd,1);1]*ones(1,nPix);  % constraint vector
Cbig   = kron(speye(nPix),C);             % matrix involved in the Hessian of the barrier
Qbig   = kron(speye(nPix),S'*S);          % matrix involved in the Hessian of the data-fitting term
StY    = S'*Y;
StS    = S'*S;                            % Hessian of the data-fitting term
minStS = min(eig(S'*S));                  % minimal eigen value of the Hessian of the data-fitting term

%%% useful functions
grad_phi = @(Z,CZ_c,mu_) -StY+StS*Z-mu_.*C'*(1./CZ_c);       % gradient of smooth term + barrier
prox     = @(Z,B,norm_B,step,prec,V) prox_metric(Z,B,...
                norm_B,step,prec,V,nRow,nCol,J,Wav,Wav_inv); % proximity operator with metric B_1

%%% initialization
X    = ones(nEnd,nPix)/(nEnd+1); % strictly feasible initial point
CX_c = C*X+c;                    % constraints should be >0
iter = 1;
time = 0;

fprintf('------------ Start PIPA -------------\n');

while time<time_max
    
    %%% store variables to plot figures
    snr_vec(iter)     = -10*log10(sum((Xtrue(:)-X(:)).^2)./sum(Xtrue(:).^2)); % signal-to-noise ratio
    obj_vec(iter)     =  obj(X);                                              % objective function 
    X_Xinf_vec(iter)  =  norm(X(:)-Xinf(:))/norm(Xinf(:));                    % normalized distance to solution
    time_vec(iter)    =  time;                                                % running time 
    if mod(iter-1,2)==0; fprintf('iter %d Xinf_vec %.2d obj %.2d time %.1f\n',iter-1, X_Xinf_vec(iter), obj_vec(iter),time); end
    
    tic        
    %%% Check stopping criterion to decrease barrier parameter    
    if iter>1 && norm(RX1(:))<eps1*mu/zeta && norm(RX2(:))<eps2*mu/zeta && RX3<eps3*mu/zeta
        mu        = mu/rho;        % decrease barrier parameter
        precision = min(1,mu*100); % precision for computing prox depends on barrier parameter
        zeta      = zeta/1.02;     % to ensure that stopping criteria decrease faster than barrier parameter
    end
    
    %%% build preconditioner
    grad_phiX     =  grad_phi(X,CX_c,mu);
    Ssmall        =  mu./CX_c.^2;
    Sbig          =  diag(sparse(Ssmall(:)));
    A             =  Qbig+Cbig'*Sbig*Cbig ;                    % preconditioner
    A_1           =  @(u) A\u ;                                % inverse of preconditioner
    N_A_1         =  1/(mu/max(max(CX_c(1:nEnd,:)))^2+minStS); % norme of inverse of preconditioner
    A_1_grad_phiX =  reshape(A_1(grad_phiX(:)),nEnd,nPix);
    
    %%% store current iterate
    Xold    = X;
    CXold_c = CX_c;
    
    %%% start backtracking
    gamma=gamma_max;
    for iter_backtrack=1:maxiter_backtrack
        if(reg>0)
            
            %%% forward-backward iteration
            [X,VV] = prox(Xold-gamma.*A_1_grad_phiX,A_1,N_A_1,gamma*reg,precision,VV);
        else
            X =  Xold-gamma.*A_1_grad_phiX;
        end
        
        %%% check if feasible point
        CX_c = C*X+c;
        if CX_c>0              

            %%% check if sufficient decrease            
            upper_bound = delta_backtrack*((X(:)-Xold(:))'*A*(X(:)-Xold(:)))/(gamma*mu);
            SXold_X     = S*(Xold-X);
            CX_CXold    = CX_c./CXold_c;
            if 0.5*sum(SXold_X(:).^2)/mu+sum(CX_CXold(:))-(nEnd+1)*nPix-sum(log(CX_CXold(:)))<=upper_bound 
                break
            end
        end
        if iter_backtrack>maxiter_backtrack; fprintf('Primal variable backstracking did not converge'); return; end
        gamma = gamma*theta_backtrack; % decrease gamma
    end

    %%% compute stopping criteria
    RX1  = Xold - X;
    RX2  = A*(Xold(:) - X(:))./gamma;
    RX3  = sum(abs(CX_c(:)./CXold_c(:)-1));
    
    time = time + toc;
    iter = iter+1;
end

%%% store variables to plot figures
snr_vec(iter)     = -10*log10(sum((Xtrue(:)-X(:)).^2)./sum(Xtrue(:).^2)); % signal-to-noise ratio
obj_vec(iter)     =  obj(X);                                              % objective function 
X_Xinf_vec(iter)  =  norm(X(:)-Xinf(:))/norm(Xinf(:));                    % normalized distance to solution
time_vec(iter)    =  time;                                                % running time 
fprintf('iter %d Xinf_vec %.2d obj %.2d time %.1f\n',iter, X_Xinf_vec(iter), obj_vec(iter),time)

%%% recap
fprintf('------------------------------------------------------\n')
fprintf('Execution time for %d iterations of PIPA: %.1f s\n', iter, time)
fprintf('Final objective function value: %.2f\n', obj_vec(end))
fprintf('Final normalized distance to solution: %.2f\n', X_Xinf_vec(end))
fprintf('Final SNR: %.2f\n', snr_vec(end))
fprintf('------------------------------------------------------\n')

function [OBJ] = obj(Z)
    ZWt               = Wav_mult(Wav,Z);
    [~,detail_coeffs] = sort_wavelet_coeffs(ZWt,nRow,nCol,nEnd,J);
    SZ                = S*Z;
    OBJ               = 0.5*norm(Y(:)-SZ(:))^2 + reg*(sum(sum(sum(abs(detail_coeffs)))));
end

end





