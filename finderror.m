function [err] = finderror(pm,pmold)
%%
%function [error] = finderror(pm,pmold)
% Computes the error in the solution. This error is computed as the L2 norm
% for all of the points in the computational domain.
%Variables
% Inputs:
%   1. pm - new solution vector ie: pm^(n+1)
%   2. pmold - old/current solution vector, ie: pm^(n)
% Outputs:
%   1. error - L2 norm of the errors in the computation domain
% Called by:
%   1. pipeFlow
% Calls:
%   None.
%Last modified on:
%   1. 7/31/14
%   2. 2015.08.30
%%
SMALL = 1e-16;                              %required to make sure that MaxErr is never infinity even if pm is zero
MaxErr = abs(pm-pmold)./(abs(pm)+SMALL);
err = max(MaxErr);