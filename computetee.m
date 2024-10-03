function [A,b,Neq] = computetee(tee,pm,rho,Nms,Neq,A,b)
%%
% function [A,b,Neq] = computetee(tee,pm,rho,Nms,Neq,A,b)
% This function will fill the A maxtrix and b vector for tees. Three
% equations are written: continuity and two momentum equations.
%Variables:
% Input:
%   1. tee: contains all of the information about a Tee. Note that this variable is a
%           structure. 
%         .D - diameter, ft.  Only constant diameter tees.
%         .in - inlet node(s), may be a 2-elelment array [#,#] for multiple
%           inlets
%         .out - outlet node(s), may be a 2-elelment array [#,#] for
%            multiple outlets
%         .run - nodes of the run (straight section)
%         .Krun - use 20 *ft(relative roughness)
%         .Kbranch - use 60 *ft(relative roughness)
%         example: tee(1)  = struct('D',dia,'in',2,'out',[3 5],'run',[2...
%            ...5],'Krun',20*ft(edval),'Kbranch',60*ft(edval));

%   2. pm: the solution vector P1...PN, mdot1...mdotN
%   3. rho: fluid desnity in slugs/ft^3
%   4. Nms: The the location of the last pressure term, so that Nms+1 is the first mdot term
%   5. Neq: The current equation that is being written
%   6. A: The maxtrix for A*pm = b;
%   7. b: the RHS vector
% Output:
%   1. A: now filled with the correct values for tee losses
%   2. b: now filled with the correct values for tee losses
%   3. Neq: Incremented twice for each minor tee in the network
%Called by:
%   1. pipeFlow
%Calls:
%   None
%Last Modified
%   2/9/16
%Change log:
%   1. 8/2/14 - first writting 2. 8/4/14 - Verified 1 inlet and 2 outlets
%   3. 8/5/14 - added 2 inlets and 1 exit - needs to be verified 4. 2/16/16
%   4. 2/9/16 -  Added ke terms.  All Ks are multiplied by the smaller
%   velocity (mass flow terms) in the energy eqations, i.e. K1*v1^2/2 and
%   K2 *v2^2/2for the head losses to the two smaller flows.
%%
for i=1:length(tee)
    C=8/(rho*pi^2);
    D=tee(i).D;
    if isscalar(tee(i).in)
        %1 input and 2 outputs
        %Now we need to figure out which K values to use for our two
        %momentum equations
        if ismember(tee(i).in,tee(i).run) && ismember(tee(i).out(1),tee(i).run)
            K1 = tee(i).Krun;
        else
            K1 = tee(i).Kbranch;
        end
        if ismember(tee(i).in,tee(i).run) && ismember(tee(i).out(2),tee(i).run)
            K2 = tee(i).Krun;
        else
            K2 = tee(i).Kbranch;
        end
        
        %Define the indicies for the pressure terms
        pin = tee(i).in;
        po1 = tee(i).out(1);
        po2 = tee(i).out(2);
        
        %Define the indicies for the masflow rate terms
        mi = tee(i).in+Nms;
        mo1 = tee(i).out(1)+Nms;
        mo2 = tee(i).out(2)+Nms;
        
        %write Linear momentun from in to out 1
        A(Neq,pin)  = 1;
        A(Neq,po1)  = -1;
        A(Neq,mi) = C*2*pm(mi)/D^4;
        A(Neq,mo1) = -C*(K1+1)*2*pm(mo1)/D^4;
        b(Neq) = C/D^4*(pm(mi)^2-(K1+1)*pm(mo1)^2);
        Neq = Neq+1;
        
        %write Linear momentun from in to out 2
        A(Neq,pin)  = 1;
        A(Neq,po2)  = -1;
        A(Neq,mi) = C*2*pm(mi)/D^4;
        A(Neq,mo2) = -C*(K2+1)*2*pm(mo2)/D^4;
        b(Neq) = C/D^4*(pm(mi)^2-(K2+1)*pm(mo2)^2);
        
        Neq = Neq+1;
        
        %Write conservation of mass
        A(Neq,mi)   = 1;
        A(Neq,mo1)  = -1;
        A(Neq,mo2)  = -1;
        b(Neq)      = 0;
        
        Neq = Neq+1;
        
    else
        %2 inputs and 1 output
        if ismember(tee(i).in(1),tee(i).run) && ismember(tee(i).out,tee(i).run)
            K1 = tee(i).Krun;
        else
            K1 = tee(i).Kbranch;
        end
        if ismember(tee(i).in(2),tee(i).run) && ismember(tee(i).out,tee(i).run)
            K2 = tee(i).Krun;
        else
            K2 = tee(i).Kbranch;
        end
        
        %Define the indicies for the pressure terms
        pin1 = tee(i).in(1);
        pin2 = tee(i).in(2);
        po   = tee(i).out;
        
        %Define the indicies for the mass flow rate terms
        mi1  = tee(i).in(1)+Nms;
        mi2  = tee(i).in(2)+Nms;
        mo   = tee(i).out+Nms;
        
        %Write the momentum equation from inlet 1 to exit
        A(Neq,pin1) =  1; %ps
        A(Neq,po)   = -1; %pe
        A(Neq,mi1) = C*(1-K1)*2*pm(mi1)/D^4;
        A(Neq,mo) = -C*2*pm(mo)/D^4;
        b(Neq) = C/D^4*((1-K1)*pm(mi1)^2-pm(mo)^2);
        
        Neq         = Neq+1;
        
        %Write the momentum equation from inlet 2 to exit
        A(Neq,pin2) =  1; %ps
        A(Neq,po)   = -1; %pe
        A(Neq,mi2)  = C*(1-K2)*2*pm(mi2)/D^4;
        A(Neq,mo)   = -C*2*pm(mo)/D^4;
        b(Neq)      = C/D^4*((1-K2)*pm(mi2)^2-pm(mo)^2);
        
        Neq         = Neq+1;
        
        %write Conservation of mass
        A(Neq,mi1)  = 1;
        A(Neq,mi2)  = 1;
        A(Neq,mo)   = -1;
        b(Neq)      = 0;
        Neq         = Neq+1;
    end
end



