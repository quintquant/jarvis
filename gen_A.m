% Function that generates a data set in the format necessary for PCA
% by using a third function whose name is given by fname to generate
% observations corresponding to different parameter values given in params_list
%
% INPUT:
%
%   fname       String giving the name of the function to be evlauated
%               to obtain the observations. The function must take
%               a single vector parameter and produce a vector or a matrix
%
%   params_list Matrix contianing parameter values. Each column gives the
%               parameter values to obtain one observation. The number
%               of columns M is therefore the number of observations.
%
% OUTPUT:
%
%   imax, jmax  Integers giving the dimensions of the output of the function.
%               The matrix A will contain each observation in a column,
%               whether the function named fname produces a vector (e.g. a
%               time sequence) or a matrix (e.g. an image). Therefore,
%               imax and jmax can be used by another function to, starting
%               from A, restore the observations to the format used by fname
%
%   A           The matrix A. Each column is an observation
%
function [A,imax,jmax]=gen_A(fname,params_list)
    M=size(params_list)(2); % Number of observations
    params=params_list(:,1);
    obs=feval(fname,params);
    [imax,jmax]=size(obs);
    A=zeros(imax*jmax,M);
    A(:,1)=reshape(obs,imax*jmax,1); % Convert image matrix to a vector format
    for i=2:M
        params=params_list(:,i);
        obs=feval(fname,params);
        A(:,i)=reshape(obs,imax*jmax,1);
    endfor
endfunction
