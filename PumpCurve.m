function [QH, dHdQ] = PumpCurve(x, y)
%
% MAE 413
%
% Establish interpolation for pump characteristic curves using the pchip
% structure.

% Calculate derivative at node points where (x,y) data are given using all the points provided.
% Establish interpolation for the derivatives
%
% Inputs are an array of Q and H values.  These arrays must be the same size.
% Each row in Q and H comprises values for one pump characteristic curve.
%
% The outputs of this function are two piecewise polynomial structures for use by PPVAL.
% The first structure is used to interplate to determine the value of H given a value of Q.
% The second structure is used to interpolation to determine the slope of the pump
%   characteristic curve.
%
% Define your pump curve(s) as part of your inputs.  Then call this
%   function:
%   [H, dHdQ] = PumpCurve(Q, H)
%
% When you need the value of pump head for a known value of
%   volumetric flow rate, call:
%   PPVAL(H, Q)
%
% When you need the slope of the pump curve at a known value of volumetric
% flow rate, call:
%   PPVAL(dHdQ, Q)
%
%NOTE: The first derivative is continuious so it does not mater which
%interval I use to determine the derivative.
%%
% Dr. Jamie Canino
% October 9th, 2014
% Based on PumpCurve.m from Dr. Brett Batson dated: October 12, 2013
%
%%
N = size(x,1);
for m = 1:N % Calculate interpolation structures for each pump curve
    QH(m,:) = pchip(x(m,:), y(m,:));
    NN = size(x,2);
    d  = zeros(NN,1);
    for i=1:NN
        if i==NN
            coef = QH(m).coefs(i-1,:);
            x0   = QH.breaks(i-1);
        else
            coef = QH(m).coefs(i,:);
            x0   = QH.breaks(i);
        end
        dcoef = polyder(coef);      %the derivative of a polynomial
        xx = x(i)-x0;
        d(i) = polyval(dcoef,xx); %evaluate the derivative at x-x0
    end
    dHdQ(m,:) = pchip(x(m,:), d);
end