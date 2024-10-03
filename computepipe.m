function [A,b,Neq] = computepipe(pipe,pm,rho,mu,gam,Nms,Neq,A,b)
%%
% function [A,b,Neq] = computepipe(elem,pm,rho,mu,Npts,Neq,A,b)
% This function will fill the A maxtrix and b vector for pipes. Two
% equations are written (continuity and momentum) for each pipe.
%Variables:
% Input:
%   1. pipe: the primary structure which contains all of the information for
% each element (pipe, bend, minor loss, etc). Note that this variable is a
% structure. See pipeFlow.m for a more detailed description.
%   2. pm: the solution vector P1...PN, mdot1...mdotN
%   3. rho: fluid desnity in slugs/ft^3
%   4. mu: fluid viscisity in lbf-s/ft^2 or a slug/(ft-s)
%   5. gam: rho*g in lbf/ft^3
%   6. j: the location of all of the pipe elemnts in the structure elem
%   7. Nms: The the location of the last pressure term, so that Nms+1 is the first mdot term
%   8. Neq: The current equation that is being written
%   9. A: The maxtrix for A*pm = b;
%   10. b: the RHS vector
% Output:
%   1. A: now filled with the correct values for pipes
%   2. b: now filled with the correct values for pipes
%   3. Neq: Incremented twice for each pipe in the network
%Called by:
%   1. pipeFlow
%Calls:
%   1.getf
%Last Modified
%   2015.08.30
%%

for i=1:length(pipe)
    %Write the momentum equation
    pin = pipe(i).start; 
    po  = pipe(i).end;
    mi  = pipe(i).start+Nms;
    mo  = pipe(i).end+Nms;

    D=pipe(i).D;
    eD=pipe(i).eD;
    mdot=pm(mo);
    f = getf(D,eD,mdot,mu);
    M = -16*f*pipe(i).L/(rho*pi^2*D^5);
    
    A(Neq,pin)   =  1; %ps
    A(Neq,po)   = -1; %pe
    A(Neq,mo)   = M*pm(mo); %mdot_e
    b(Neq)      = gam*pipe(i).dz + M/2*pm(mo)^2;
    Neq         = Neq+1;
    
    %write continuity
    A(Neq, mi)  = 1;
    A(Neq, mo)  = -1;
    b(Neq)      = 0;
    Neq         = Neq+1;
end
