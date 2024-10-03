function [A,b,Neq] = computepump(pump,pm,Nms,Neq,A,b)

for i=1:length(pump)
    %write the momentum equation
    pin = pump(i).start;
    po = pump(i).end;
    mi  = pump(i).start+Nms;
    %mo  = pump(i).end+Nms;

    [QH, dHdQ] = PumpCurve(x, y);
    dPlocal = ppval(QH, x(i));
    dPdmdot  = ppval(dHdQ, x(i));

    %Write the momentum equation
    A(Neq,pin)  =  1; %pin
    A(Neq,po)   = -1; %po
    A(Neq,mi)   = dPdmdot*pm(mi); %mdot i
    %A(Neq,mo)   = ; %mdot o
    b(Neq)      = dPlocal-(dPdmdot*pm(mi));
    Neq         = Neq+1;


    %Write the continuity equation
    A(Neq, mi)  = 1;
    A(Neq, mo)  = -1;
    b(Neq)      = 0;
    Neq         = Neq+1;
end