% ========================================================================================================== %
% This software was written in 2020 through 2022                                                             %
% Kirk H. Bevan <kirk.bevan@mcgill.ca> & Botong Maio <botong.miao@mail.mcgill.ca>                            %
%                                                                                                            %
% LICENCE                                                                                                    %
% Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License.                    %
% This software is distributed without any warranty.                                                         %
% You should have received a copy of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0           %
% International Public License along with this software.                                                     %
% If not, see <https://creativecommons.org/licenses/by-nc-sa/4.0/>.                                          %
%                                                                                                            %
% DISCLAIMER                                                                                                 %
% The authors and publishers make no warranties about the software, and disclaim liability for all uses of   %
% the software, to the fullest extent permitted by applicable law.  The authors do not recommend use of this %
% software for any purpose.  It is made freely available, solely to clarify points made in their analysis of %
% photoelectrochemical (PEC) devices.  When using or citing the software, you should not imply endorsement   %
% by the authors.                                                                                            %
%                                                                                                            %
% DESCRIPTION                                                                                                %
% Numerically simulates the current-voltage properties at semiconductor-liquid interface based on Gartner's  %
% model and Peter's model. In this code, we assume n-type tantalum nitride, a high supporting electrolyte    %
% concentration, and an akaline liquid region. Hence, changes to system parameters for alternate materials   %
% and electrolytes, beyond those values provided, may require direct changes to the code to obtain           %
% physically sound solutions.                                                                                %
%                                                                                                            %
%                                                                                                            %
% References                                                                                                 %
%   1. L.M. Peter, "Energetics and kinetics of light-driven oxygen evolution at semiconductor electrodes:    %
%      the example of hematite" J, Solid State Electrochem. 17, 315 (2013).                                  %
%   2. L.M. Peter and K.G.U. Wijayantha, "Photoelectrochemical water splitting at semiconductor electrodes:  % 
%      Fundamental problems and new perspectives" ChemPhysChem, 15, 1983 (2014).                             %
%   3. W.W. Gartner, "Depletion-layer photoeffects in semiconductors" Phys. Rev., 116, 84 (1959). 	     %
%   4. Kirk H. Bevan, Botong Miao, and  Asif Iqbal, "SLJCompact: A Semiconductor-Liquid Junction Solver for  % 
%      Rapid Band Diagram Insights into Photoelectrochemical Devices" (2022).                                %
%                                                                                                            %
% ========================================================================================================== %

close all;
clear all;

% ===========================================================================================================
% ==================....Parameters, General Variables, and Operator Definitions Start....====================

% -----------------------------------------------------------------------------------------------------------
% Input the fundamental parameters
q = 1.602e-19;                           % Charge on an electron [coul] 
kB = 1.38e-23;                           % Boltzmann constant [J/K]
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m] 
T = 300;                                 % Temperature in Kelvin 
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]
I0 = 6.2063e16*(100^2);                  % Incident photon flux [m^-2 s^-1]
Vt = kB*T/q;                             % Thermal voltage

% -----------------------------------------------------------------------------------------------------------
% Calculate & set semiconductor parameters
% Some benchmarks here... http://www.ioffe.ru/SVA/NSM/Semicond/
% For tantalum nitride see... https://doi.org/10.1007/s40243-016-0083-z
epsr = 17;                                                      % Relative dielectric constant, see DOI: 10.1039/c4cp05616g
eps = epsr*eps0;                                                % Overall dielectric constant [F/m]
ktr = 10;                                                       % Hole transfer rate [s^-1]
r = 2.5e-9;                                                     % Hole capture radius, see R. F. Pierret "Advanced Semiconductor Fundamentals" (2002)
sigmap = pi*r^2;                                                % Capture cross section [m^2]
vth = 6.24e4;                                                   % Thermal velocity [m/s], see R. F. Pierret "Advanced Semiconductor Fundamentals" (2002)
alpha = 3.4e7;                                                  % Absorption coefficient 400 nm [m^-1], see https://doi.org/10.1007/s40243-016-0083-z
Eg = 2.2;                                                       % Band gap [eV], https://doi.org/10.1007/s40243-016-0083-z
ND = 1e23;                                                      % n-type: Donor density [#/m^3]
Lp = sqrt((2*eps*Eg)/(q*ND));                                   % Approximate hole diffusion length [m]
Wp = 2.5e-9;                                                    % Hole capture depth [m]
Nc = 1.1106e+26;                                                % Effective density of conduction band states
Nv = 1.6815e+26;                                                % Effective density of valence band states
p0 = (Nc*Nv*exp(-Eg/Vt))/ND;                                    % Bulk hole concentration
Vbi = 0.25;                                                     % Built in voltage
carrier_max = 5.00e+28;                                         % Density of atoms in silicon

% -----------------------------------------------------------------------------------------------------------
% Setup the external grid
dx = Wp/100;
x_min = -dx*round(3*Lp/dx);
x = [x_min:dx:0];

% -----------------------------------------------------------------------------------------------------------
% Apply external bias
% This model assumes that the reverse bias is positive, and assumes the flatband potential as the zero point. 
dV = Vt/4;                     % Potential Bias Step
Vmax = 2.3;                    % Maximum potential under illumination
Vmin = 0;                      % Minimum potential under illumination
Vbias=[Vmin:dV:Vmax];          % Potential range under illumination

% ====================....Parameters, General Variables, and Operator Definitions End....====================
% ===========================================================================================================


% ===========================================================================================================
% ============================....Peter's Model & Gartner's Model Start....==================================
tic

% -----------------------------------------------------------------------------------------------------------
% Calculate photocurrent based on Peter's model

W0=sqrt(2*eps/(q*ND));                                % Depletion layer thickness for a potential of 1 volt across it
Wd=W0.*sqrt(Vbias);                                  % Depletion layer width 
g=I0.*(1-exp(-(alpha.*Wd))./(1+alpha*Lp));           % Hole flux 
phi_sc=q.*ND.*Wd.^2./(2*eps);                        % Potential drop in space charge region
n_surf=ND.*exp(-(Vbias./Vt));                        % Surface electron concentration 
krec=sigmap.*vth.*n_surf;                             % Surface recombination rate 
p_surf=g./(krec+ktr);                                 % Surface hole concentration 
J_rec=q.*krec.*p_surf;                                % Recombination current density 
J_photo_P=q.*g-J_rec;                                 % Peter's Photocurrent 

% -----------------------------------------------------------------------------------------------------------
% Calculate carrier concentrations and current-voltage properties via modified form of L. Peter's model

Wd_bi = sqrt((2*eps*Vbi)/(q*ND));
phi_b = ((q*ND)/(2*eps))*(x+Wd_bi).^2;
[val,Wd_bi_index] = min(abs(x+Wd_bi));
phi_b(1:Wd_bi_index) = 0;
p_eq = p0*exp(phi_b/Vt);
for ii=1:length(Vbias)
  Vbb = Vbias(ii);
  phi_b = ((q*ND)/(2*eps))*(x+Wd(ii)).^2;
  [val,Wd_index] = min(abs(x+Wd(ii)));
  phi_b(1:Wd_index) = 0;
  p_s(ii) = (p_surf(ii) + sum(p_eq.*exp(-abs(x)/Wp))*dx) ...
            /sum(exp((phi_b-phi_b(end))/Vt).*exp(-abs(x)/Wp)*dx);
  n_s(ii) = ND*exp(-Vbb/Vt);
  p_light = p_s(ii)*exp((phi_b-phi_b(end))/Vt);
  p_dark = p0*exp(phi_b/Vt);
  J_dark(ii) = q*ktr*sum((p_dark-p_eq).*exp(-abs(x)/Wp))*dx;
  J_light(ii) = q*ktr*sum((p_light-p_eq).*exp(-abs(x)/Wp))*dx;
  if(p_s(ii) > 0.001*p_dark(end))
    J_light(ii) = J_light(ii) + J_dark(ii);   % Assume two competing statistics when the quasi-Fermi levels are close to closing
  end
  DeltaV_F(ii) = -Vt*[log(p_dark(end)/Nv)-log(p_s(ii)/Nv)];
  p_light_list{ii} = p_light;
  p_dark_list{ii} = p_dark;
  phi_b_list{ii} = phi_b;
  n_list{ii} = ND*exp(-phi_b/Vt);
  Ec = phi_b - phi_b(end);
  Ev = Ec-Eg;
  light_Fn{ii} = Ec(1)*ones(1,length(x)) + Vt*log(ND/Nc);
  %light_Fp{ii} = Ev - Vt*log(p_light/Nv);   % Flat throughout
  light_Fp{ii} = light_Fn{ii} - DeltaV_F(ii).*[1-1./(1+exp((x+Wd(ii))/Wp))];   % Merged approximation
  light_Ec{ii} = Ec;
  light_Ev{ii} = Ev;
  light_ps(ii) = sum(p_light.*exp(-abs(x)/Wp))*dx;
  light_ns(ii) = sum(n_list{ii}.*exp(-abs(x)/Wp))*dx;
end

Vonset = -Vt*log(ktr)+Vt*log(pi*(r^2)*vth*ND);


% -----------------------------------------------------------------------------------------------------------
% Calculate the quasi-Fermi level comparatively through the current-voltage results
for ii=1:length(Vbias)
   Vbb_dark(ii) = interp1(J_dark,Vbias,J_light(ii));
   dV_currents(ii) = Vbb_dark(ii) - Vbias(ii);
end

toc
% =============================....Peter's Model & Gartner's Model End....===================================
% ===========================================================================================================


% -----------------------------------------------------------------------------------------------------------
% Current-voltage figures
figure;
set(gcf,'Position',[100 200 400 900])
subplot(4,1,1);
plot(Vbias,J_dark/(q*I0),'b'); ylim([-0.1 2]);  hold on;
plot(Vbias,J_light/(q*I0),'r',Vbias,J_photo_P/(q*I0),'k--'); ylim([0 2]);
xlabel('V (V_{Ref}-V_{FB})');
ylabel('J_{p,ht} (Arb.)');
xlim([0 Vmax]);
grid on;

subplot(4,1,2);
semilogy(Vbias,light_ps/(Wp*carrier_max),'r',Vbias,light_ns/(Wp*ND),'k');
xlabel('V (V_{Ref}-V_{FB})');
ylabel('Conc. (Arb.)');
legend('surface holes','surface electrons','location','SouthEast');
xlim([0 Vmax]);
ylim([1e-7 2]);

subplot(4,1,3);
semilogy(Vbias,J_rec./(q*g),'r',Vbias,(q*g)./(q*g),'k');
xlabel('V (V_{Ref}-V_{FB})');
ylabel('Flux (Arb.)');
legend('R_{eh,flux}/G_{flux}','G_{flux}/G_{flux}','location','SouthEast');
xlim([0 Vmax]);
ylim([1e-5 2]);

subplot(4,1,4)
plot(Vbias,DeltaV_F,'k',Vbias,dV_currents,'g'); 
xlabel('V (V_{Ref}-V_{FB})');
ylabel('\DeltaV_{F} (Arb.)');
xlim([0 Vmax]);
ylim([0.8 1.7]);
legend('Analytical','Light-Dark Current Offset');

% -----------------------------------------------------------------------------------------------------------
% Band diagram plotting at Vonset
figure;
[val,ip] = min(abs(Vbias-Vonset));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x/1e-9,light_Ec{ip},'r'); hold on;
plot(x/1e-9,light_Ev{ip},'r'); hold on;
plot(x/1e-9,light_Fn{ip},'r--'); hold on;
plot(x/1e-9,light_Fp{ip},'r--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(light_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Illumination Band Diagram | %1.2f V',Vbias(ip));
title(s);
