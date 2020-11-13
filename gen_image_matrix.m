% Function to convert an image (e.g. a neutron scattering image)
% back to matrix form when it is given as a single vector:
%
function matrix=gen_image_matrix(vector,imax,jmax)
    matrix=reshape(vector,imax,jmax);
endfunction
