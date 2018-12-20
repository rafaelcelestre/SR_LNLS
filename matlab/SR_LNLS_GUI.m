%% SR - calculations (v1.0.0b)

% Campinas, May 26th 2015
% X-Ray Optics Group - Brazilian Light Source (LNLS)
% Author: Rafael Celestre
% rafael.celestre@lnls.br
% + 55 19 3521 1285

% This program is part of the bachelor thesis entitled 'Código para computação de 
% parâmetros de dispositivos de inserção com aplicações às simulações de ray-tracing 
% para linhas de luz de raios-X duros', written by Rafael Celestre (2014).

% This is a beta version, please report to the author for bugs and missbehaves.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DON'T EDIT BELLOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = SR_LNLS_GUI(varargin)
% SR_LNLS_GUI MATLAB code for SR_LNLS_GUI.fig
%      SR_LNLS_GUI, by itself, creates a new SR_LNLS_GUI or raises the existing
%      singleton*.
%
%      H = SR_LNLS_GUI returns the handle to a new SR_LNLS_GUI or the handle to
%      the existing singleton*.
%
%      SR_LNLS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SR_LNLS_GUI.M with the given input arguments.
%
%      SR_LNLS_GUI('Property','Value',...) creates a new SR_LNLS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SR_LNLS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SR_LNLS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SR_LNLS_GUI

% Last Modified by GUIDE v2.5 09-Oct-2014 14:01:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SR_LNLS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SR_LNLS_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SR_LNLS_GUI is made visible.
function SR_LNLS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SR_LNLS_GUI (see VARARGIN)

% Choose default command line output for SR_LNLS_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clc
disp('SR - calculations (v1.0.0b)');


disp('Campinas, May 26th 2015');
disp('X-Ray Optics Group - Brazilian Light Source (LNLS)');
disp('Authors: Rafael Celestre & Bernd Meyer');
disp('rafael.celestre@lnls.br');
disp('+ 55 19 3521 1285');
%% Pre-loading SIRIUS parameters
% Storage Ring Parameters
set(handles.r_energy,'String','3.00');
set(handles.r_current,'String','500');
set(handles.gamma,'String','5871');
set(handles.emittance,'String','0.272e-9');
set(handles.coupling,'String','1.00');
set(handles.es,'String','0.083');
set(handles.beta_x,'String','1.5');
set(handles.beta_y,'String','1.4');
% Electron Beam Parameters
set(handles.hor_siz,'String','20.0988');
set(handles.hor_div,'String','13.3992');
set(handles.ver_siz,'String','1.94173');
set(handles.ver_div,'String','1.38695');
% Undulator Parameters
set(handles.B,'String','1.15');
set(handles.lambda_u,'String','19.0');
set(handles.gap,'String','5');
set(handles.N,'String','105');
set(handles.num,'String','-');
set(handles.L,'String','1.995');
set(handles.K,'String','2.04');
set(handles.k_min,'String','0.4');
set(handles.h1_energy,'String','1.461');
set(handles.Pt,'String','7.512')
% Simulation parameters
set(handles.distance,'String','-')
set(handles.Ei,'String','-')
set(handles.Ef,'String','-')
set(handles.step,'String','-')

% UIWAIT makes SR_LNLS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SR_LNLS_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Global Variables

global E_r i_r nat_emittance coupling es beta_x beta_y B lambda_u N num distance Pt
global L k k_min gamma E K lambda gap norma delta_eh delta_ev div_eh div_ev nat_size   
global nat_div delta_und div_und delta_tv div_tv delta_th div_th ph_flux brilliance
global coh_flux coh_frac P_char E_char F_flux F_brilliance F_cohflux Mag S Z Bs Bz 
global S_hd Z_hd Bs_hd Bz_hd X Xs PS_x Xp Y Ys Yp PS_y P h1 Ei Ef step

%% Directory

format long
curr_folder=pwd;
addpath(curr_folder);

%% Set up %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage Ring Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_r = str2double(get(handles.r_energy,'String'));
i_r = str2double(get(handles.r_current,'String'))/1000;
nat_emittance = str2double(get(handles.emittance,'String'));
coupling = str2double(get(handles.coupling,'String'));
es = str2double(get(handles.es,'String'))/100;
beta_x = str2double(get(handles.beta_x,'String'));
beta_y = str2double(get(handles.beta_y,'String'));
%%%%%%%%%%%%%%%%%%%%%%%%
% Undulator Parameters %
%%%%%%%%%%%%%%%%%%%%%%%%
B = str2double(get(handles.B,'String'));
lambda_u = str2double(get(handles.lambda_u,'String'))/1000;
gap = str2double(get(handles.gap,'String'));
N = str2double(get(handles.N,'String'));
num = str2double(get(handles.num,'String'));
k_min = str2double(get(handles.k_min,'String'));
distance = str2double(get(handles.distance,'String'));
%%%%%%%%%%%%%%%%
% Energy Range %
%%%%%%%%%%%%%%%%
Ei   = str2double(get(handles.Ei,'String'));
Ef   = str2double(get(handles.Ef,'String'));
step = str2double(get(handles.step,'String'));
%%%%%%%%%%%%%%%%
% Calculations %
%%%%%%%%%%%%%%%%
P_char = get(handles.photon_char,'Value');
E_char = get(handles.elec_char,'Value');
F_flux = get(handles.ph_flux,'Value');
F_brilliance = get(handles.brilliance,'Value');
F_cohflux = get(handles.coh_flux,'Value');
Mag = get(handles.magnetic,'Value');
%%%%%%%%%%%%%%%%
% Output files %
%%%%%%%%%%%%%%%%
Graph = get(handles.graphs,'Value');
Harmonics = get(handles.harmonic_txt,'Value');
Calculations = get(handles.calculation_txt,'Value');

%% Inicialization %%

%%%%%%%%%%%%%%%%%%%%%%
% Physical Constants %
%%%%%%%%%%%%%%%%%%%%%%
e  = 1.60217653E-19;                % Electron charge [C]
m  = 9.10938188E-31;                % Electron mass [kg]
c  = 2.99792458E+08;                % Light speed [m.s^-1]
h  = 6.62606957E-34;                % Planck Constant [m^2.kg.s^-1]
mu = 4*pi*1E-7;                     % Magnetic constant [N.A^-2]
epsilon = 1/(mu*c^2);               % Electric constant [F.m^-1]
%%%%%%%%%%%%%%%%%%%%%%
% Photon beam Energy %
%%%%%%%%%%%%%%%%%%%%%%
gamma = (E_r*e*1E9)/(m*c^2);        % Lorentz Factor
step = step/1000;
E = Ei:step:Ef;                     % Energy [keV]
lambda = h*c./(E*e*1E3);            % X-Ray wavelenght [m] 
%%%%%%%%%%%%%%%%%%%%%
% E-beam parameters %
%%%%%%%%%%%%%%%%%%%%%
coup = coupling/100;
emittance_y = nat_emittance * coup/(coup+1);   
emittance_x = nat_emittance * 1/(1+coup);        
% gamma functions
gamma_x = 1/beta_x;
gamma_y = 1/beta_y;
% electron beam size [µm]
delta_eh  = sqrt(beta_x*emittance_x)*1E6;
delta_ev  = sqrt(beta_y*emittance_y)*1E6;
% electron beam divergence [µrad]
div_eh  = sqrt(gamma_x*emittance_x)*1E6;
div_ev  = sqrt(gamma_y*emittance_y)*1E6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ondulator Auxiliary parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L =lambda_u*N;                      % Undulator lenght [m]
k = e*B*lambda_u/(2*pi*m*c);        % Max deflexion parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First harmonic energy [keV] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = 0.0095*E_r^2./(lambda_u*(1+(k^2)/2));
%%%%%%%%%%%%%%%%%%%%
% Total power [kW] %
%%%%%%%%%%%%%%%%%%%%
Z0 = c*mu;                          % Free space impedance
Pt = 1E-3*N*Z0*i_r*e*2*pi*c*(gamma*k)^2 /(6*lambda_u);

%% Display basic calculations %%
enable_elec = isnan(delta_eh*div_eh*delta_ev*div_ev);
enable_undu = isnan(gamma*L*k*h1);

if isnan(gamma)== 0
    set(handles.gamma,'String',gamma);
end
if  enable_elec == 0
    set(handles.hor_siz,'String',delta_eh);
    set(handles.hor_div,'String',div_eh);
    set(handles.ver_siz,'String',delta_ev);
    set(handles.ver_div,'String',div_ev);
end
if enable_undu == 0
    set(handles.L,'String',L);
    set(handles.K,'String',k);
    set(handles.h1_energy,'String',h1);
end
if enable_undu*i_r == 0
	set(handles.Pt,'String',Pt);
end

%% Calculations %%

enable = 0;
enable_sr = isnan(E_r*i_r*nat_emittance*coup*es*beta_x*beta_y);
enable_un = isnan(B*lambda_u*gap*N*num*k_min);
enable_sm = isnan(Ei*Ef*step*distance);
if P_char||E_char||F_flux||F_brilliance||F_cohflux||Mag == 1
    enable_cl = 0;
else
    enable_cl = 1;
end
enable = enable_sr+enable_un+enable_sm+enable_cl;
if enable == 0
    path = uigetdir();                  % allows to choose which directory to work at
    curr_folder=pwd;                    % saves the current adress
    cd(path);                           % loads the chosen path
    if F_cohflux == 1                   
        F_brilliance = 1;        
    end
    if F_brilliance == 1                
        P_char = 1;
        F_flux = 1;
    end
    %%%%%%%%%%%%%%%%%%%
    % Magnetic fields %
    %%%%%%%%%%%%%%%%%%%
    if Mag == 1
        l_u = lambda_u*1000;                % lambda_u [mm]
        Bo = B*cosh(pi*gap/l_u);            % B field at the magnet
        ku = 2*pi/l_u;                      % spatial wavenumber
        s = 2*l_u;                          % path twice as long as lambda_u
        % 3D plot-friendly
        [S,Z]= meshgrid(-s/2:(s/170):s/2,-gap/2:(gap/30):gap/2);
        Bz = Bo*cosh(ku*Z).*cos(ku*S)./cosh(pi*gap/l_u); 
        Bs = Bo*sinh(ku*Z).*sin(ku*S)./cosh(pi*gap/l_u);
        % 2D (HD matrix)
        [S_hd,Z_hd]= meshgrid(-s/2:0.01:s/2,-gap/2:0.01:gap/2);
        Bz_hd = Bo*cosh(ku*Z_hd).*cos(ku*S_hd)./cosh(pi*gap/l_u);
        Bs_hd = -Bo*sinh(ku*Z_hd).*sin(ku*S_hd)./cosh(pi*gap/l_u);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Electron beam Phase-space %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if E_char == 1
         if delta_eh>div_eh             % chooses the longest axis
            x_range  = delta_eh*3;      % sets the axis to be 3x sigma
            xp_range = x_range;         % makes axis symmetric
        else
            x_range  = div_eh*3;
            xp_range = x_range;
        end
        if delta_ev>div_ev              % chooses the longest axis
            y_range  = delta_ev*3;      % sets the axis to be 3x sigma
            yp_range = y_range;         % makes axis symmetric
        else
            y_range  = div_ev*3;
            yp_range = y_range;
        end
        x_step   = 5000;                % HD image
        xp_step  = x_step;
        y_step   = xp_step;
        yp_step  = y_step;
        % Phase space calcualtions
        x  = 1/(delta_eh*sqrt(2*pi)); 
        xp = 1/(div_ev*sqrt(2*pi));
        stepx  = x_range/x_step; 
        stepxp = xp_range/xp_step;
        y  = 1/(delta_ev*sqrt(2*pi));
        yp = 1/(div_eh*sqrt(2*pi));
        stepy  = y_range/y_step; 
        stepyp = yp_range/yp_step;
        [Xs, Xp] = meshgrid(-x_range/2:stepx:x_range/2, -xp_range/2:stepxp:xp_range/2);
        [Ys, Yp] = meshgrid(-y_range/2:stepy:y_range/2, -yp_range/2:stepyp:yp_range/2);
        % Horizontal phase-space
        PS_x = 1000*i_r*x*xp*exp(-2*(Xs/delta_eh).^2).*exp(-2*(Xp/div_eh).^2);
        % Vertical phase-space
        PS_y = 1000*i_r*y*yp*exp(-2*(Ys/delta_ev).^2).*exp(-2*(Yp/div_ev).^2);
        % E-beam profile
        [X, Y] = meshgrid(-x_range/2:stepx:x_range/2, -x_range/2:stepx:x_range/2);
        P = 1000*i_r*x*y*exp(-2*(X/delta_eh).^2).*exp(-2*(Y/delta_ev).^2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Photon beam calculations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    if P_char||F_flux||F_brilliance||F_cohflux == 1
        K	 = zeros(num,size(E,2));
        norma = zeros(num,size(E,2));
        %%%%%%%%%%%%
        % K matrix %
        %%%%%%%%%%%%
        for i = 1:1:num
            n = 2*(i-1)+1;                  % Odd harmonics
            K(i,:) = sqrt((4*n*lambda*gamma^2)/lambda_u - 2);       
            for j = 1:size(E,2)
                Re = real(K(i,j));
                Im = imag(K(i,j));
                if Im~=0
                    K(i,j)=0;
                elseif Re>k||Re<k_min       % Upper and lower limits for K
                    K(i,j)=0;
                end
                if K(i,j)~=0    
                    norma(i,j)=1;            % Auxiliar matrix for selecting data
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%
        % Characteristics %
        %%%%%%%%%%%%%%%%%%%
        if P_char == 1
            P_char = get(handles.photon_char,'Value');
            nat_size  = zeros(num,size(E,2));
            nat_div   = zeros(num,size(E,2));
            delta_und = zeros(num,size(E,2));
            div_und   = zeros(num,size(E,2));
            delta_tv  = zeros(num,size(E,2));
            div_tv    = zeros(num,size(E,2));
            delta_th  = zeros(num,size(E,2));
            div_th    = zeros(num,size(E,2));
            for i = 1:1:num                
                n = 2*(i-1)+1;                  % Odd harmonics
                epsilon_e = 2*pi*n*N*es;         
                % Undulator Natural Size
                nat_size(i,:) = (1E6)*sqrt(2*L*lambda)/(4*pi); % [µm]
                % Undulator Natural Divergence
                nat_div(i,:) = (1E6)*sqrt(lambda/(2*L)); %[µm]
                % Undulator source size
                delta_und(i,:) = Qs(epsilon_e) * nat_size(i,:); % [µm]
                % Undulator source divergence
                div_und(i,:) = Qa(epsilon_e) * nat_div(i,:); % [µrad]
                % Vertical photon size
                delta_tv(i,:) = sqrt((delta_ev^2) + delta_und(i,:).^2); % [µm]
                % Vertical photon divergence
                div_tv(i,:) = sqrt((div_ev^2) + div_und(i,:).^2); % [µrad]
                % Horizontal photon size
                delta_th(i,:) = sqrt((delta_eh^2) + delta_und(i,:).^2); % [µm]
                % Horizontal photon divergence
                div_th(i,:) = sqrt((div_eh^2) + div_und(i,:).^2);  % [µrad]        
            end
            % Setting the lower and upper Energy limits
            nat_size    = nat_size .*norma;
            nat_div     = nat_div  .*norma;
            delta_tv    = delta_tv .*norma;
            div_tv      = div_tv   .*norma;
            delta_th    = delta_th .*norma;
            div_th      = div_th   .*norma;
            % Beam Propagation     
            if distance ~= 0
                delta_tv    = 1E6*distance*2*tan((div_tv*1E-6)/2) + delta_tv;
                delta_th    = 1E6*distance*2*tan((div_th*1E-6)/2) + delta_th;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Photon Flux [Ph/s/1% bandwidth]  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if F_flux == 1
            F_flux = get(handles.ph_flux,'Value');
            Bessel  = zeros(num,size(E,2));
            arg     = zeros(num,size(E,2));
            Fn      = zeros(num,size(E,2));
            Qn      = zeros(num,size(E,2));
            ph_flux = zeros(num,size(E,2));
            bw = 0.1/100;                       % Bandwidth
            alfa = (e^2)/(2*c*h*epsilon);       % fine structure constant 
            konst = 0.5*alfa*pi*bw/e;
            for i=1:num
                n = 2*(i-1)+1;                  % Odd harmonics
                arg(i,:) = n*(K(i,:).^2)./(4 + 2 * K(i,:).^2);
                Bessel(i,:)= besselj((n-1)/2, arg(i,:)) - besselj((n+1)/2, arg(i,:));
                Fn(i,:) = (n*K(i,:).*Bessel(i,:)./(1 + 0.5 * K(i,:).^2)).^2;
                Qn(i,:) = Fn(i,:).*(1 + 0.5 * K(i,:).^2)/n;       
                ph_flux(i,:) = konst*N*Qn(i,:)*i_r;    
            end
            ph_flux = ph_flux.*norma;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Brilliance [Ph/s/mm²/mrad²/0.1%BW %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if F_brilliance == 1
            F_brilliance = get(handles.brilliance,'Value');
            brilliance = ph_flux./(1E-12*(4*pi^2)*div_tv.*div_th.*delta_tv.*delta_th);
            for j=1:size(E,2)
                for i=1:num
                    if norma(i,j)==0
                        brilliance(i,j)=0;
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Coherence calculation %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        if F_cohflux == 1
            F_cohflux = get(handles.coh_flux,'Value');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Coherent flux [Ph/s/1% bandwidth] %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:num
                coh_flux(i,:) = 0.25*1E6*brilliance(i,:).*((lambda*1E3).^2);
            end
            coh_flux = coh_flux.*norma;
            %%%%%%%%%%%%%%%%%%%%%
            % Coherent fraction %
            %%%%%%%%%%%%%%%%%%%%%
            for i=1:num
                coh_frac(i,:) = 1E8*((lambda*1E3/(4*pi)).^2)./(1E-12*div_tv(i,:)...
                    .*div_th(i,:).*delta_tv(i,:).*delta_th(i,:));
            end
            coh_frac = coh_frac.*norma;
        end
    end
    f_txt()
end
%% Saving generated data
if Graph == 1
	plot_graph()
end
if Harmonics == 1
	txt_harm()
end
if Calculations == 1
    txt_calc()
end

%% Return to code folder %%
cd(curr_folder);                        % loads the original folder address
uiwait(msgbox('Calculation finished',' ','modal'));

% --- Executes during object creation, after setting all properties.
function r_energy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r_current_Callback(hObject, eventdata, handles)
% hObject    handle to r_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_current as text
%        str2double(get(hObject,'String')) returns contents of r_current as a double


% --- Executes during object creation, after setting all properties.
function r_current_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function emittance_Callback(hObject, eventdata, handles)
% hObject    handle to emittance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of emittance as text
%        str2double(get(hObject,'String')) returns contents of emittance as a double


% --- Executes during object creation, after setting all properties.
function emittance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to emittance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coupling_Callback(hObject, eventdata, handles)
% hObject    handle to coupling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coupling as text
%        str2double(get(hObject,'String')) returns contents of coupling as a double


% --- Executes during object creation, after setting all properties.
function coupling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coupling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_x_Callback(hObject, eventdata, handles)
% hObject    handle to beta_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_x as text
%        str2double(get(hObject,'String')) returns contents of beta_x as a double


% --- Executes during object creation, after setting all properties.
function beta_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_y_Callback(hObject, eventdata, handles)
% hObject    handle to beta_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_y as text
%        str2double(get(hObject,'String')) returns contents of beta_y as a double


% --- Executes during object creation, after setting all properties.
function beta_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function B_Callback(hObject, eventdata, handles)
% hObject    handle to B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B as text
%        str2double(get(hObject,'String')) returns contents of B as a double


% --- Executes during object creation, after setting all properties.
function B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_u_Callback(hObject, eventdata, handles)
% hObject    handle to lambda_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_u as text
%        str2double(get(hObject,'String')) returns contents of lambda_u as a double


% --- Executes during object creation, after setting all properties.
function lambda_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gap_Callback(hObject, eventdata, handles)
% hObject    handle to gap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gap as text
%        str2double(get(hObject,'String')) returns contents of gap as a double


% --- Executes during object creation, after setting all properties.
function gap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in photon_char.
function photon_char_Callback(hObject, eventdata, handles)
% hObject    handle to photon_char (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of photon_char


% --- Executes on button press in elec_char.
function elec_char_Callback(hObject, eventdata, handles)
% hObject    handle to elec_char (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of elec_char


% --- Executes on button press in ph_flux.
function ph_flux_Callback(hObject, eventdata, handles)
% hObject    handle to ph_flux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ph_flux


% --- Executes on button press in brilliance.
function brilliance_Callback(hObject, eventdata, handles)
% hObject    handle to brilliance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of brilliance


% --- Executes on button press in coh_flux.
function coh_flux_Callback(hObject, eventdata, handles)
% hObject    handle to coh_flux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of coh_flux


% --- Executes on button press in magnetic.
function magnetic_Callback(hObject, eventdata, handles)
% hObject    handle to magnetic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of magnetic


% --- Executes on button press in graphs.
function graphs_Callback(hObject, eventdata, handles)
% hObject    handle to graphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of graphs


% --- Executes on button press in harmonic_txt.
function harmonic_txt_Callback(hObject, eventdata, handles)
% hObject    handle to harmonic_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of harmonic_txt


% --- Executes on button press in calculation_txt.
function calculation_txt_Callback(hObject, eventdata, handles)
% hObject    handle to calculation_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calculation_txt



function num_Callback(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num as text
%        str2double(get(hObject,'String')) returns contents of num as a double


% --- Executes during object creation, after setting all properties.
function num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function distance_Callback(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of distance as text
%        str2double(get(hObject,'String')) returns contents of distance as a double


% --- Executes during object creation, after setting all properties.
function distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ei_Callback(hObject, eventdata, handles)
% hObject    handle to Ei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ei as text
%        str2double(get(hObject,'String')) returns contents of Ei as a double


% --- Executes during object creation, after setting all properties.
function Ei_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ei (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ef_Callback(hObject, eventdata, handles)
% hObject    handle to Ef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ef as text
%        str2double(get(hObject,'String')) returns contents of Ef as a double


% --- Executes during object creation, after setting all properties.
function Ef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_Callback(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step as text
%        str2double(get(hObject,'String')) returns contents of step as a double


% --- Executes during object creation, after setting all properties.
function step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function es_Callback(hObject, eventdata, handles)
% hObject    handle to es (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of es as text
%        str2double(get(hObject,'String')) returns contents of es as a double


% --- Executes during object creation, after setting all properties.
function es_CreateFcn(hObject, eventdata, handles)
% hObject    handle to es (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k_min_Callback(hObject, eventdata, handles)
% hObject    handle to k_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_min as text
%        str2double(get(hObject,'String')) returns contents of k_min as a double


% --- Executes during object creation, after setting all properties.
function k_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r_energy_Callback(hObject, eventdata, handles)
% hObject    handle to r_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_energy as text
%        str2double(get(hObject,'String')) returns contents of r_energy as a double
