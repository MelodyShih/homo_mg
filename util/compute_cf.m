% Antti Hannukainen / 4.3.2018 / Otaniemi
%
%-----------------------------------------------
%
% Fit a line to logarithm of the error using LSQ-method.
% The slope of this line is called "contraction factor".
%
% Datapoints with relative error smaller than 1e-9 are omitted.
%


function [cf,cof] = compute_cf(error)

M = max(error);
b = log(error(error/M > 1e-9));

Niter = 1:length(b);
A = [Niter(:) ones(size(Niter(:))) ];
cof = A\b(:);
cf = cof(1);
