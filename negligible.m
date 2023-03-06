% Adapted from ncomp_tools -- see below. Free license given on this software:

function result = negligible(X,Y,tol)

% Return true if X is "negligible" with respect to Y with tolerance 'tol'
% We must be careful with zero elements of Y...

%************************************************************************%
%                                                                        %
% ncomp_tools: L. Barnett, C. L. Buckley and S. Bullock (2009)           %
%                                                                        %
% [*] See http://www.secse.net/ncomp/ncomp1.pdf                          %
%                                                                        %
% IF YOU USE THIS CODE PLEASE CITE THE ABOVE AS:                         %
%                                                                        %
%   L. Barnett, C. L. Buckley and S. Bullock (2009)                      %
%   On Neural Complexity and Structural Connectivity                     %
%   Physical Review E (in review)                                        %
%                                                                        %
%************************************************************************%

% J Lizier - At some point octave was complaining about this line, but it seems ok now
result = (nnz(abs(X) > tol*eps(max(abs(Y),eps))) == 0);

% An alternative (not as nice but workaround) was :
%result = (nnz(abs(X) > tol*eps) == 0);
