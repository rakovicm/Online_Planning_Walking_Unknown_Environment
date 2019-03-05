function s = vmod(a)
%
% scalar = vmod( [3x3]_vector )
% returns the module of a spatial vector
% or a vector of modules of some spacial vecotrs vecor

s = sqrt(diag(a'*a));
