% Function to project a given set of observations onto a given (incomplete)
% set of given principal compoments:
%
function s=scores(Xo,U,R)
   s=U(:,1:R)'*Xo;
endfunction
