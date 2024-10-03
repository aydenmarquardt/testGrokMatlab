%%
% pipeFlow.m - main program - HW 3 Fixed
% Jamie Canino, Ph.D.
% modified by J. Koch
% Last Updated:
% 8/20/23

% This code will compute the pressure, in psf, and flow rate, in slugs/s,
% at nodes in a piping network.
% This code utilizing a structure array which contains all of the
% information required for constructing the piping network.


%Calls:
% 1. ft.m - computes the fully-turbulent friction factor
% 2. computepipe.m - constructs the A and b matricies for all of the pipes
% 3. computetee.m - constructs the A and b matricies for all of the tees
% 4. boundaryConditions.m - constructs the A and b matricies for the
% pressure boundary conditions.
% 5. finderror.m - computes the error for convergence.

%Called by:
% none

%Variables:
% 1. pipe - structure with information for 
%   .L = length in ft
%   .D = dimater in ft
%   .start = the starting location of the pipe
%   .end = the ending location of the pipe
%   .eD = the relative roughness of the pipe, unitless
%   .dz = Zend-Zstart in ft
% 2. tee - structure with information for 
%   .D = dimater in ft
%   .in = the inlet(s) of the tee
%   .out = the outlet(s) of the tee
%   .run = the two tee points that are across from eachother
%   .Krun = the C value for a run, usually 20
%   .Kbranch = the C value for the branch, usually, 60
% 3. bc - is a structure array for the boundary conditions
%   .type - the type of boundary condition. Right now only 'pressure' is a
%   valid option
%   .val - the constant value for the boundary condition (in psf for type =
%   pressure)
%   .loc - the point at which the boundary condition should be applied.
% 4. pm - vecture that contains first the pressures and then the massflow 
%         rates at each point
% 5. rho - fluid density in slugs/ft^3
% 6. mu - fluid viscosity in lbf-s/ft^2 or a slug/(ft-s)
% 7. g - gravitational constant, ft/s^2
% 8. gam - specific weight of the fluid, lbf/ft^3
% 9. N - the number of equation required
% 10. Nms - the location of the last pressure term in the vector pm
% 11. tol - the tolerance for convergence
% 12. MAXiter = the maximum number of times the code should iterate
% 13. A - the matrix that must be inverted to solve for pressure and 
%           mass flow rates, as in A*pm = b. 
% 14. b - the right hand side vector

%Reference: Batson, B., Pipe_flow_rev_21.m, September 2013
%%
clear all
close all
clc
format shortg
%%
% This section of the code is the inputs
%dia = 4.026/12;
%edval = 4.47e-4;
%pipe(1) = struct('L',50,'D',dia,'start',1,'end',2,'eD',edval,'dz',0);
%bc(1)   = struct('type','pressure','val',(2+14.7)*144,'loc',1);
%bc(2)   = struct('type','pressure','val',14.7*144,'loc',2);

%dia = 0.936/12;
% dia=0.957/12;
% edval = 2e-4/dia;
% pipe(1) = struct('L',6,'D',dia,'start',1,'end',2,'eD',edval,'dz',1);
% bc(1)   = struct('type','pressure','val',(16)*144,'loc',1);
% bc(2)   = struct('type','mfr','val',0.07639,'loc',1);

%minor(1) = struct('Di',4/12,'Do',2/12,'start',1,'end',2,'Ki',0);
%bc(1) = struct('type','pressure','val',(14.7)*144,'loc',1);
%bc(2) = struct('type','pressure','val',(8)*144,'loc',2);

dia1 = .23/12; %small D in ft
dia2 = .35/12; %big D in ft
edval1 = 0.00006/dia1; %relative roughness for pipes of D1 and D2
edval2 = 0.00006/dia2;

mdot = [2.1 164.3 270.9 366.9 554.7 663.7 757.3 849.1 928 1009.1 1094.4]; %in gpm
dP = [108.7 109.1 108.3 106.6 101.9 97.2 92.5 86.5 81 74.2 66.5]; %in feet water column

x = mdot*0.004324; %provides x input in slugs/s for PumpCurve.m
y = dP*5.20233; %provides y input in psf for PumpCurve.m

%pump(1) = struct() %necessary in structure to call computepump.m
pipe(1) = struct('L',2.5/12,'D',dia1,'start',1,'end',2,'eD',edval1,'dz',0);
tee(1) = struct('D',dia1,'in',2,'out',[3,4],'run',[3,4],'Krun',20*ft(edval1),'Kbranch',30*ft(edval1));
pipe(2) = struct('L',11.5/12,'D',dia1,'start',3,'end',5,'eD',edval1,'dz',11.5/12);
minor(1) = struct('Di',dia1,'Do',dia2,'start',5,'end',6,'Ki',30*ft(edval1)+(1-(dia1/dia2)^2)^2); %this models the elbow + contraction as 1 element
pipe(3) = struct('L',26.75/12,'D',dia2,'start',6,'end',7,'eD',edval2,'dz',0);
minor(2) = struct('Di',dia2,'Do',dia2,'start',7,'end',8,'Ki',30*ft(edval2));
pipe(4) = struct('L',3/12,'D',dia2,'start',8,'end',18,'eD',edval2,'dz',-3/12);
minor(6) = struct('Di',dia2,'Do',dia1,'start',18,'end',19,'Ki',.5*(1-(dia1/dia2)^2)*sqrt(sin(90)));
pipe(9) = struct('L',2.5/12,'D',dia1,'start',19,'end',9,'eD',edval1,'dz',-2.5/12);
tee(2) = struct('D',dia1,'in',[9,10],'out',11,'run',[9,10],'Krun',20*ft(edval1),'Kbranch',30*ft(edval1));
pipe(6) = struct('L',22.75/12,'D',dia1,'start',14,'end',10,'eD',edval1,'dz',22.75/12);
minor(4) = struct('Di',dia1,'Do',dia1,'start',15,'end',14,'Ki',30*ft(edval1));
pipe(7) = struct('L',28.25/12,'D',dia1,'start',16,'end',15,'eD',edval1,'dz',0);
minor(5)  = struct('Di',dia1,'Do',dia1,'start',17,'end',16,'Ki',30*ft(edval1));
pipe(8) = struct('L',16.75/12,'D',dia1,'start',4,'end',17,'eD',edval1,'dz',16.75/12);
pipe(5) = struct('L',9.5/12,'D',dia1,'start',11,'end',12,'eD',edval1,'dz',0);
minor(3) = struct('Di',dia1,'Do',dia1,'start',12,'end',13,'Ki',30*ft(edval1));

bc(1) = struct('type','pressure','val',22.17*144,'loc',1);
bc(2) = struct('type','pressure','val',17.8*144,'loc',13);

%%
% Constants for the code
Relax = 1;
rho = 1.94;         %slugs/ft^3, at 60 F (Hodge and Taylor, Table B-2)
mu  = 2.34e-5;  %lbf-s/ft^2 or a slug/(ft-s), at 60 F (Hodge, Table B-2)
g   = 32.174;       %ft/s^2
tol = 1e-8;         %tolerance
MAXiter = 40;

%%
%The code requires that all structures be defined. In order to clean up the
%inputs (and prevent errors) if a certain element is not in the piping
%network (say a minor loss) then the corresponding structure (minor in
%this case) is set to an empty array.
if exist('pipe','var')==0
    pipe = [];
end
if exist('tee','var')==0
    tee = [];
end
if exist('minor','var')==0
    minor = [];
end
%%
% Initalize some variables
gam = rho*g;
N = 2*length(pipe)+3*length(tee)+length(bc)+2*length(minor);     %Total number of equations
%% error trap
if mod(N,2)~=0
    error('Why isn"t the number of equations even?');
end
%%
Nms = N/2;          % location of the last pressure term, so that Nms+1 is the first mdot term
A = zeros(N,N);
b = zeros(N,1);
pm = 0.1*ones(N,1);
j = find(strcmpi({bc.type},'pressure')==1);
pm(1:N/2) = max([bc(j).val]);               %Set all of the pressures to the inlet pressure
%%

%write some stuff to the command line
NN = length(pipe)+length(tee);
fprintf('Number of piping elements = %3.0f,  A is a %d x %d matrix \n', NN, N, N);
fprintf('\n');
%%
pmold = pm;
err   = 1;      %Initialzie the error to make the loop start
iter  = 0;
while err>=tol && iter<MAXiter
    iter = iter +1;
    Neq = 1;                             %Neq will give us the row we are currently working on, ie the equation number
    if ~isempty(pipe)
        [A,b,Neq] = computepipe(pipe,pm,rho,mu,gam,Nms,Neq,A,b);
    end
    if ~isempty(tee)
        [A,b,Neq] = computetee(tee,pm,rho,Nms,Neq,A,b);
    end
    if ~isempty(minor)
        [A,b,Neq] = computeminor(minor,pm,rho,Nms,Neq,A,b);
    end
    [A,b,Neq]=boundaryConditions(bc,Neq,A,b,Nms);
    [A,b,Neq] = computepump(pump,pm,Nms,Neq,A,b);
    pm = A\b;
    pm = Relax*pm + (1-Relax)*pmold;
    err = finderror(pm,pmold);
    error_out(iter) = err;
    pmold=pm;
    fprintf('After %d iterations, MaxErr = %8.6e \n', iter, err);
end
%%
%Write the solution
fprintf('\n');
fprintf('   Node     mdot            Q         p          rho            mu \n');
fprintf('          [slugs/s]       [gpm]     [psia]     [slug/cf]     [lbf*s/ft^2] \n');
for i=1:Nms
    fprintf(' %4.0f  %10.5f  %13.5f   %7.5f     %5.4f       %6.5e \n', ...
        i, pm(Nms+i),pm(Nms+i)/rho/0.13368*60,pm(i)/144, rho, mu);
end
fprintf('\n');

%% troubleshooting outputs
%pm
%semilogy(1:iter,error_out,'s-')









