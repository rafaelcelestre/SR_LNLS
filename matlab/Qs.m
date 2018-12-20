function [ y ] = Qs( x )

% Function used to introduce energy spread effects to the ondulator photon beam.
% For more info, refer to B. C. Meyer, "Nota técnica - OPT 010-13 - Photon Size
% and Divergence," Campinas, 2013.

x_s=x/4;
y = 2*(Qa(x_s)^(2/3));

end

