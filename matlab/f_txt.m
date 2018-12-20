function [ ] = f_txt( )

% Function to write log file. No input is required.

global E_r i_r nat_emittance coupling es beta_x beta_y B lambda_u N gamma
global L k k_min delta_eh delta_ev div_eh div_ev Pt h1 Ei Ef step distance


DateVector = clock;                  % gets the current date
DateString = datestr(DateVector);    % transforms date to string
fid = fopen('parameters.txt', 'w');     % opens a txt file to be written
fprintf(fid,'%s \n \n', DateString);
fprintf(fid,'Planar Undulator Calculations\n\n');
% Storage ring parameters
fprintf(fid,'Storage Ring Parameters:\n');
fprintf(fid,'\tRing Energy:         %.2f [GeV]\n', E_r);
fprintf(fid,'\tRing Current:        %d [mA]\n', i_r*1E3);
fprintf(fid,'\tGamma:               %.4f\n', gamma);
fprintf(fid,'\tNatural emittance:   %d\t[m.rad]\n', nat_emittance);
fprintf(fid,'\tCoupling factor:     %d\t\t\t\t[%%]\n', coupling);
fprintf(fid,'\tEnergy spread:       %d\t[%%]\n', es*100);
fprintf(fid,'\tHor_beta:            %d\t[m]\n', beta_x);
fprintf(fid,'\tVer_beta:            %d\t[m]\n\n', beta_y);
% Electron beam parameters
fprintf(fid,'E-beam parameters (sigma):\n');
fprintf(fid,'\tHor_size:            %d\t[um]\n',  delta_eh);
fprintf(fid,'\tHor_div:             %d\t[urad]\n', div_eh);
fprintf(fid,'\tVer_size:            %d\t[um]\n', delta_ev);
fprintf(fid,'\tVer_div:             %d\t[urad]\n\n', div_ev);
% Undulator parameters
fprintf(fid,'Undulator parameters:\n  ');
fprintf(fid,'\tMagnetic field:      %d\t[T]\n', B);
fprintf(fid,'\tMagnetic period:     %d\t\t\t\t[mm]\n', 1000*lambda_u);
fprintf(fid,'\tNumber of periods:   %d\t\t\t[ ]\n', N);
fprintf(fid,'\tUndulator Lenght:    %d\t[m]\n', L);
fprintf(fid,'\t1st harmonic energy: %d\t[keV]\n', h1);    
fprintf(fid,'\tK max:               %d\t[ ]\n', k);
fprintf(fid,'\tK min:               %d\t[ ]\n', k_min);
fprintf(fid,'\tTotal Power:         %d\t[kW]\n\n', Pt);
% Simulation parameters
fprintf(fid,'Simulation parameters:\n  ');
fprintf(fid,'\tDist. to the source: %.2f\t[m]\n', distance);
fprintf(fid,'\tInitial energy:      %.2f\t[keV]\n', Ei);
fprintf(fid,'\tFinal Energy:        %.2f\t[keV]\n', Ef);
fprintf(fid,'\tEnergy step:         %.2f\t[ev]\n', step*1000);
fclose(fid);

end