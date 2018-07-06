
%====================================================================
% Algorithm 12 on page 161 from
% Pustelnik, Nelly. Méthodes proximales pour la résolution de problèmes 
% inverses: application à la tomographie par émission de positrons. 
% Diss. Université Paris-Est, 2010.
%====================================================================



function NormL = CalculNorme( L, x )
% These Nelly : Algo 12 page 161

% Calcule ||   L || 

% % Initialisation
% 
% p = 1 + 10^(-6) ;
% pnew = 1 ;
% 
% n = 1 ;
% 
% epsilon = 10^(-8) ;
% 
% nmax = 20;
% 
% cond = abs( pnew-p ) / pnew 

% Initialisation

pnew = 1 ;

n = 1 ;

epsilon = 0.01 ;

nmax = 20;

cond = 1; 

% Iterations

while ( cond >= epsilon && n < nmax)
    
    xnew = L(x) ;
    
    p = pnew ;
    pnew = norm(xnew) / norm(x) ;
    
    cond = abs(  pnew-p ) / pnew ;
    
    x = xnew;
    n = n+1 ;
end

% Calcul de la norme

 NormL = p;

end


