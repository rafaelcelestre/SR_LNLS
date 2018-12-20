function [] = txt_harm()
% txt_harm aims to generate a .txt file for each simulated harmonic. files are to be 
% separeted by a ' '. from left to right the included data is: energy[keV], lambda[A]
% K[], sigma_hor[um], sigma_ver[um], sigma_hor_p[urad], sigma_ver_p[urad],
% ph_flux[ph/s/0.1bw], brilliance[ph/s/mm^2/mrad^2/0.1bw]

global num E lambda K delta_tv div_tv delta_th div_th ph_flux brilliance coh_flux norma
global P_char F_flux  F_brilliance F_cohflux

pchar = 0;
fflux = 0;
fbrilliance = 0;
fcohflux = 0;

if P_char == 0
    delta_tv = ones(num,size(E,2)).*norma;
    div_tv   = ones(num,size(E,2)).*norma;
    delta_th = ones(num,size(E,2)).*norma;
    div_th   = ones(num,size(E,2)).*norma;
    pchar    = 1;
end
    
if F_flux == 0
    ph_flux     = ones(num,size(E,2)).*norma;
    brilliance  = ones(num,size(E,2)).*norma;
    coh_flux    = ones(num,size(E,2)).*norma;
    fflux       = 1;
end

if F_brilliance == 0
    brilliance  = ones(num,size(E,2)).*norma;
    fbrilliance = 1;
end

if F_cohflux == 0
    coh_flux    = ones(num,size(E,2)).*norma;
    fcohflux    = 1;
end

for i=1:1:num
    E_aux = E(1,:)';
    lambda_aux = (lambda(1,:)')*1E10;
    E_aux = E_aux.*(norma(i,:)');
    E_aux(E_aux==0)=[];
    lambda_aux = lambda_aux.*(norma(i,:)');
    lambda_aux(lambda_aux==0)=[];
    n     = 2*(i-1)+1;
    K_a   = K(i,:)';
    K_a(K_a==0)=[];
    d_tv  = delta_tv(i,:)';
    d_tv(d_tv==0)=[];
    di_tv = div_tv(i,:)';
    di_tv(di_tv==0)=[];
    d_th  = delta_th(i,:)';
    d_th(d_th==0)=[];
    di_th = div_th(i,:)';
    di_th(di_th==0)=[];
    if pchar == 1
        d_tv = d_tv*NaN;
        di_tv = di_tv*NaN;
        d_th = d_th*NaN;
        di_th = di_tv*NaN;
    end
    p_f   = ph_flux(i,:)';
    p_f(p_f==0)=[];
    if fflux == 1
        p_f = p_f*NaN;
    end
    b     = brilliance(i,:)';
    b(b==0)=[];
    if fbrilliance == 1
        b = b*NaN;
    end
    c_f   = coh_flux(i,:)';
    c_f(c_f==0)=[];
    if fcohflux == 1
        c_f = c_f*NaN;
    end
    data_out = [E_aux,lambda_aux,K_a,d_tv,di_tv,d_th,di_th,p_f,b,c_f];
    filename = sprintf('harm_%d.txt', n);
    fid = fopen(filename,'wt');
    fprintf(fid, '%s %s %s %s %s %s %s %s %s %s\n', 'energy[keV]','lambda[A]','K',...
        'sig_v[um]','div_v[urad]','sig_h[um]','div_h[urad]','ph_flux[ph/s/0.1%bw]',...
        'brilliance[ph/s/mm2/mrad2/0.1%bw]','coh.flux[ph/s/0.1%bw]');
    dlmwrite(filename, data_out,'delimiter', ' ', '-append')
    fclose(fid);
end

end

