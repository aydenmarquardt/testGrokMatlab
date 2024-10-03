function [A,b,Neq]=boundaryConditions(bc,Neq,A,b,Nms)
%%
%function [A,b,Neq]=boundaryConditions(bc,Neq,A,b)
% Function fills A and b for the boundary conditions. Currently this code
% is only written for a constant static pressure boundary condition.
%Variables
% Inputs:
%   1. bc - structure containing the information for the boundary condition
%   2. Neq - the current equation that is being written
%   3. A - A*pm = b
%   4. b - the right-hand-side vector
% Outputs:
%   1. A - now filled
%   2. b - now filled
%   3. Neq - incremented for each boundary condition that was written
% Called by:
%   1. pipeFlow
% Calls:
%   None.
%Last modified on:
%   1. 7/31/14
%   2. 2015.08.30
%   3. 2023.08.14
%%
for i = 1: length(bc)
    switch bc(i).type
        case 'pressure'
            A(Neq,bc(i).loc) = 1;
            b(Neq)           = bc(i).val;
            Neq = Neq+1;
        case 'mfr' %mass flow rate (slugs/s)
            A(Neq,bc(i).loc+Nms) = 1;
            b(Neq)               = bc(i).val;
            Neq = Neq+1;
        otherwise
            error('\n  Unrecognized bc type in bc(%0.0f).',i)
    end
end



