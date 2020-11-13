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
    X=A-ones(size(A)(1),1)*mean(A); % MATLAB: X=A-ones(size(A,1),1)’*mean(A);
    % Singular value decomposition:
    [V,lambda,junk] = svd(X'*X);
    % Principal components:
    U = X*V*lambda^(-.5);
    % Principal component scores:
    S = U'*X;
    % Data for scree plot:
    scree_data=diag(lambda);
endfunction
