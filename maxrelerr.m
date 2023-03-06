% Adapted from ncomp_tools -- see below. Free license given on this software:

function mre = maxrelerr(X,Y)

% Return maximum relative error of X with respect to Y
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

mre = max(max(abs(X)./max(abs(Y),eps)));
