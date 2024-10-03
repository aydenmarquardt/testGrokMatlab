function [A,b,Neq] = computeminor(minor,pm,rho,Nms,Neq,A,b)
%this function should calculate the head loss for a pipe network part that
%contains minor losses
%Variables:
% Input
%   1. minor: A structured variable that contains the beginning and ending
%  diameter
%   2. pm: the solution vector P1...PN, mdot1...mdotN
%   3. rho: fluid desnity in slugs/ft^3
%   4. Nms: The the location of the last pressure term, so that Nms+1 is the first mdot term
%   5. Neq: The current equation that is being written
%   6. A: The maxtrix for A*pm = b;
%   7. b: the RHS vector
% Output
%   1. A: now filled with the correct values for minor loss elements
%   2. b: now filled with the correct values for minor loss elements
%   3. Neq: Incremented twice for each pipe in the network
% Called By
%   1. pipeFlow
% Calls
%   1. ft
% Last Modified
%   2024.9.11
%%

for i=1:length(minor)
    %write the momentum equation
    pin = minor(i).start;
    po = minor(i).end;
    mi  = minor(i).start+Nms;
    mo  = minor(i).end+Nms;
    
    Di = minor(i).Di; %diameter in
    Do = minor(i).Do; %diameter out
    Ki = minor(i).Ki; %minor loss constant, taken at inlet

    Mi = 16/(rho*pi^2*Di^4); %this is just a constant, being precalculated for speed NOTE: NOT MASS FLOW RATE
    Mo = 16/(rho*pi^2*Do^4);

    %Write the momentum equation
    A(Neq,pin)  =  1; %pin
    A(Neq,po)   = -1; %po
    A(Neq,mi)   = Mi*pm(mi)*(1-Ki); %mdot i
    A(Neq,mo)   = -Mo*pm(mo); %mdot o
    b(Neq)      = (Mi/2) * pm(mi)^2 - (Mo/2) * pm(mo)^2 - (Ki*Mi/2) * pm(mi)^2;
    Neq         = Neq+1;


    %Write the continuity equation
    A(Neq, mi)  = 1;
    A(Neq, mo)  = -1;
    b(Neq)      = 0;
    Neq         = Neq+1;
end

