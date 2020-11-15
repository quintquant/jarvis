% Reconstruction of X using only the 1st, 2nd, ..., and Rth principal components:
function Xr=recons(U,S,R)
    M=size(U)(2);
    for m=1:M
        Xr(:,m)=U(:,1:R)*S(1:R,m); 
    endfor
endfunction
