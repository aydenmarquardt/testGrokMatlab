%attempts to derive the friction factor
%mdot = mass flow rate (slugs/s or kg/s)
%d = pipe diameter (ft or m)
%e = relative roughness (unitless)
%mu = dynamic viscosity (lb*s/ft2 or N*s/m2
%usage: ffinder(mdot,d,e,mu)
%ensure all units match!!!
%either lbf, slug, ft or N, kg, m

function f_guess = getf(d,e,mdot,mu)

Re = mdot/(mu*d*pi/4); %this was hard for me to derive for some reason
%if you're unsure if you're doing the units right, make sure the Re comes out unitless

if Re<2300 %if it's fully laminar, take the easy way out
    f_guess = 64/Re;
    return
else
    f_guess = (-1.8*log10((e/3.7)^1.11 + (6.91/Re)))^(-2); %use the explicit (inaccurate) method as a starting guess
    %often the explicit method is good enough and the program only needs
    %to loop once
    loop = 1;
    error = 1;
    f_new = zeros(10,1); %preallocate memory for the array of past iterations, which is kept for debugging
    while error>.0001
    
        f_new(loop) = (-2*log10((e/3.7)+(2.51/(Re*sqrt(f_guess)))))^(-2); %here's the implicit method for getting f

        error = abs(f_new(loop)-f_guess)/f_guess; %take the difference between iterations as a percentage of the previous iteration's value
    
        f_guess = f_new(loop); %if the error is still too much, this will set us up for a second iteration. If the error is acceptable, this is our output term
        
        loop = loop+1;
        if loop>=10
            return %panic stop if it goes too long (no convergence)
        end
    
    end
end


end