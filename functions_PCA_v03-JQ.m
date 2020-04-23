1==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differences between v03 and v02:
%
% Added function scores(...) to determine the decomposition of a given 
% image on given principal components (WORK IN PROGRESS)
%
% Differences between v01 and v02:
%
% In the definition of function PCA(A), we have replaced the line 
%    X=A-mean(A); 
% with 
%    X=A-ones(size(A)(1),1)*mean(A);
% The former syntax works in Octave 4.4.1 (installed on workstation) but is
% not backward-compatible with Octave 3.4.3 (installed on tor), where 
% it returns an error. The latter works with both versions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
%%%%%%%%

% Function to carry out Principal Component Analysis of given data
%
% INPUT:
%
%   data    N x M matrix with containing M observations of an object consisting 
%           of N data points in each case. data(n,m) is the nth data point 
%           for the mth observation. 
%
%           Examples: 
%
%           1) An observation could be an N-pixel two-dimensional image, such
%           as a neutron-scattering pattern, and 
%           each data point could be a pixel in that image; different
%           observations could correpsond to measurements of the 
%           neutron scattering pattern in the same material under differnet 
%           experimental conditions, 
%           e.g. temperature or applied magnetic field, or even for different 
%           materials within the same family (M sets of experimental conditions 
%           or materials in total). The neutron-scattering
%           images obtained in different observations would have to
%           correspond to the same vlaues of the scattering vector (qx,qy).
%
%           2) An observation 
%           could be a time sequence of a real variable, e.g. moun asymmetry
%           vs time, and each data point 
%           could be the value of the variable at a particular time, out
%           of N times given in a feixed sequence. 
%           Each new observation could correspond to obtaining the same
%           time sequence under differnet experimental conditions, e.g.
%           at a different temperature, with a total of M values of temperature.
%           The time stamps t would
%           have to be the same for different observations.
%
% OUTPUT:
%
%   X       N x M matrix dientical to A but with the background substracted. 
%           More precisely, the average of all values for each 
%           observation has been substracted from that observation.
%
%   U       N x M Matrix containing the M principal components. The mth column  
%           of U is the mth principal component. Note that the total number
%           of principal components is the same as the number of observations.
%           The principal components are given in order of increasing variance
%           in the data set that they capture: U(:,1) is the PC capturing the
%           greatest amount of variance in the data set, U(:,2) is the PC that 
%           is second best at capturring the variance in the data set, etc.
%
%   S       M x M (square) matrix containing the scores of the principal 
%           components: score(i,j) is the score of the jth observation 
%           on the ith principal component, i.e. the amplitude corresponding 
%           to the ith PC in the decomposition of the jth observation 
%           in principal components:
%
%           X(:,m) = U*S(:,m) 
%
%           i.e. X(:,m) = U(:,1)*S(1,m) + U(:,2)*S(2,m) + ... + U(:,M)*S(M,m)
%
%           Because the principal components, i.e. the columns of U,
%           come given in order of increasing variance, is is possible to
%           get a good reconstruction of X using a smaller number R < M of 
%           principal components:
%
%           X(:,m) approx.  = U(:,1:R)*S(1:R,m)
%
%                           = U(:,1)*S(1,m) + U(:,2)*S(2,m) + U(:,R)*S(R,m)
%   

function [X,U,S,scree_data]=PCA(A)
    % Substract the mean from each observation:
    %   Example: if each observation is a function of time A(t)
    %   and different observations were taken at differnet tempeatures,
    %   here we are substracting the time-average of A(t) for a given
    %   temperature from the values of A(t) obtained at different times for
    %   that temperature.
    X=A-ones(size(A)(1),1)*mean(A); 
    % Singular value decomposition:
    [V,lambda,junk] = svd(X'*X);
    % Principal components:
    U = X*V*lambda^(-.5);   % Eq 13 in Robbie Twyman's report
    % Principal component scores:
    S = U'*X;          
    % Data for scree plot:
    scree_data=diag(lambda);
endfunction

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

% TEST:
%
%octave:> N=      2;
%octave:> Gamma=  0.4;
%octave:> Delta=  0.0;
%octave:> h=      [0.00,0.50,1.00,1.50];
%octave:> T=      0.02;
%octave:> q_max=  4.5;
%octave:> Nk=6;
%octave:> lengths=zeros(7,1);
%octave:> lengths(1)=length(N);
%octave:> lengths(2)=length(Gamma);
%octave:> lengths(3)=length(Delta);
%octave:> lengths(4)=length(h);
%octave:> lengths(5)=length(T);
%octave:> lengths(6)=length(q_max);
%octave:> lengths(7)=length(Nk);
%octave:> params_list=zeros(7,max(lengths));
%octave:> params_list(1,:)=N;
%octave:> params_list(2,:)=Gamma;
%octave:> params_list(3,:)=Delta;
%octave:> params_list(4,:)=h;
%octave:> params_list(5,:)=T;
%octave:> params_list(6,:)=q_max;
%octave:> params_list(7,:)=Nk;
%octave:> [A,imax,jmax]=gen_A("gen_pic_mat_wrapped",params_list)

% Function to convert an image back to matrix form when it is given 
% as a single vector:
%
function matrix=gen_image_matrix(vector,imax,jmax)
    matrix=reshape(vector,imax,jmax);
endfunction

% Reconstruction of X using only the 1st, 2nd, ..., and Rth principal components:
function Xr=recons(U,S,R)
    M=size(U)(2);
    for m=1:M
        Xr(:,m)=U(:,1:R)*S(1:R,m)
    endfor
endfunction

% Function to project a given set of observations onto a given (incomplete) set of given principal compoments (WORK IN PROGRESS)
%
function s=scores()
%   WORK IN PROGRESS
endfunction
