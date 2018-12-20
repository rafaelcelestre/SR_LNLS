function [] = txt_calc()
% txt_calc aims to generate a .txt file for each calculation performed. Files are to be 
% separeted by a ' '.


global E lambda delta_tv div_tv delta_th div_th ph_flux brilliance coh_flux
global P_char F_flux  F_brilliance F_cohflux

%%%%%%%%%%%%%%%%%%%
% Characteristics %
%%%%%%%%%%%%%%%%%%%
if P_char == 1
    data_out = [E',lambda',delta_tv'];
    filename = sprintf('vertical size.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','sig_v[um]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
    
    data_out = [E',lambda',div_tv'];
    filename = sprintf('vertical divergence.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','div_v[urad]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
    
    data_out = [E',lambda',delta_th'];
    filename = sprintf('horizontal size.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','sig_h[um]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
    
    data_out = [E',lambda',div_th'];
    filename = sprintf('horizontal divergence.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','div_h[urad]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Photon Flux [Ph/s/1% bandwidth]  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if F_flux == 1
    data_out = [E',lambda',ph_flux'];
    filename = sprintf('photon flux.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','ph_flux[ph/s/0.1%bw]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brilliance [Ph/s/mm²/mrad²/0.1%BW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if F_brilliance == 1
    data_out = [E',lambda',brilliance'];
    filename = sprintf('brilliance.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','brilliance[ph/s/mm2/mrad2/0.1%bw]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coherent flux [Ph/s/1% bandwidth] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if F_cohflux == 1
    data_out = [E',lambda',coh_flux'];
    filename = sprintf('coherent flux.txt');
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s\n', 'energy[keV]','lambda[A]','coh.flux[ph/s/0.1%bw]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
end