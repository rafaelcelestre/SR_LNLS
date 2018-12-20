function [ y ] = Qa( x )

% Function used to introduce energy spread effects to the ondulator photon beam.
% For more info, refer to B. C. Meyer, "Nota técnica - OPT 010-13 - Photon Size
% and Divergence," Campinas, 2013.
if x == 0
    y = 1;
else
    y = sqrt(2*(x^2)/(-1+exp(-2 * x^2)+sqrt(2*pi)*x*erf(sqrt(2)*x)));
end
if y<=1
    y = 1;
end
end

