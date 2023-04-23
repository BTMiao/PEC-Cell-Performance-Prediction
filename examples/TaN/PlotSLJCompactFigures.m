% ========================================================================================================== %
% This software was written in 2022 by                                                                       %
% Kirk H. Bevan <kirk.bevan@mcgill.ca>, Botong Miao, and Asif Iqbal ("the authors").                         %
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
% o Example script for plotting the output from SLJCompactModel.m pertaining to the primary reference below. %
%                                                                                                            %
% Primary Reference                                                                                          %
%   1. Kirk H. Bevan, Botong Miao, and  Asif Iqbal, "SLJCompact: A Semiconductor-Liquid Junction Solver for  %
%      Rapid Band Diagram Insights into Photoelectrochemical Devices" (2022).                                %
%                                                                                                            %
% Closely Related References                                                                                 %
%    2.  K. Li, B. Miao, W. Fa, R. Chen, J. Jin, K. H. Bevan and D. Wang, “Evolution of Surface Oxidation on %
%        Ta3N5 as Probed by a Photoelectrochemical Method”, ACS Appl. Mater. Interfaces 13, 17420 (2021).    %
%    3.  B. Miao, K. Sangare, A. Iqbal, B. Marsan, and K. H. Bevan, “Interpreting interfacial                %
%        semiconductor-liquid capacitive characteristics impacted by surface states: a theoretical and       %
%        experimental study of CuGaS2”, Phys. Chem. Chem. Phys, 22, 19631 (2020).                            %
%    4.  B. Miao, A. Iqbal, and K. H. Bevan, "Utilizing Band Diagrams to Interpret the Photovoltage and      %
%        Photocurrent in Photoanodes: A Semiclassical Device Modeling Study", J. Phys. Chem. C 123, 28593    %
%        (2019).                                                                                             %
%    5.  A. Iqbal, S. Yuan, Z. Wang, and K. H. Bevan "The Impact of Bulk Trapping Phenomena on the Maximum   %
%        Attainable Photovoltage of Semiconductor-Liquid Interfaces”, J. Phys. Chem. C 122, 23878 (2018).    %
%    6.  A. Iqbal and K. H. Bevan “The impact of boundary conditions on calculated photovoltages and         %
%        photocurrents at photocatalytic interfaces”, MRS Communications 8, 466 (2018).                      %
%    7.  A. Iqbal and K. H. Bevan, "On Simultaneously Solving the Photovoltage and Photocurrent at           %
%        Semiconductor-Liquid Interfaces”, J. Phys. Chem. C 122, 30 (2018).                                  %
%    8.  A. Iqbal, Md. S. Hossain and K. H. Bevan, "The role of relative rate constants in determining       %
%        surface state phenomena at semiconductor–liquid interfaces", Phys. Chem. Chem. Phys. 18, 29466      %
%        (2016).                                                                                             %
%                                                                                                            %
% ========================================================================================================== %

close all;
clear all;
beep off;


% ========================================================================================================== %
% Fundamental Constants
q = 1.602e-19;                           % Charge on an electron [coul]
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m]
kB = 1.38e-23;                           % Boltzmann constant [J/K]
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]

% ==========================================================================================================
% Inputs
light_filename = 'Light_kp1_tau-8.mat';   % Change these files to plot for alternate runs
dark_filename = 'Dark_kp1.mat';
Bias = 0.605;
logflag = 1;

% ==========================================================================================================
% Load calculated data
load(light_filename);
Vt = sys.T*kB/q;
light_sys = sys;
light_data = jv_data
load(dark_filename);
dark_sys = sys;
dark_data = jv_data;
light_data.Jp = light_data.Jp/q;   % Normalize
dark_data.Jp = dark_data.Jp/q;     % Normalize

% ==========================================================================================================
% Post processing & selection calculations

% ------------------------------------------------------------------
% Equilibrium Band Bending 
Vbb_eq = abs(light_data.phi_b_eq(end)-light_data.phi_b_eq(1))
x = light_data.x;
[val,Nsolid] = min(abs(x))
x_solid = x(1:Nsolid);

% ------------------------------------------------------------------
% Compute the quasi-Fermi levels and splitting at each illuminated bias
light_Fp = {};
light_Fn = {};
light_Ec = {};
light_Ev = {};
light_dF = [];
ps = [];
ns = [];
for ii=1:length(light_data.Vbb)
  Ec = light_data.phi_b{ii}(1:Nsolid) - light_data.phi_b{ii}(end);
  Ev = Ec-sc.Eg;
  light_Fn{ii} = Ec(1)*ones(1,Nsolid) + Vt*log(sc.ND/sc.Nc);
  light_Fp{ii} = Ev - Vt*log(light_data.p{ii}/sc.Nv);
  light_Ec{ii} = Ec;
  light_Ev{ii} = Ev;
  light_dF(ii) = light_Fn{ii}(Nsolid)-light_Fp{ii}(Nsolid);   % Check dFn splitting at large biases
  p_temp = (light_data.p{ii}).*(exp(-abs(x_solid)/inter.Wp));
  light_ps(ii) = sum((p_temp(1:end-1)+p_temp(2:end))*.5.*diff(x_solid));
  n_temp = (light_data.n{ii}).*(exp(-abs(x_solid)/inter.Wp));
  light_ns(ii) = sum((n_temp(1:end-1)+n_temp(2:end))*.5.*diff(x_solid));
  %light_ns(ii) = light_data.n{ii}(Nsolid);
end

% ------------------------------------------------------------------
% Compute the dark bias band diagram properties
dark_Ef = {};
dark_Ec = {};
dark_Ev = {};
for ii=1:length(light_data.Vbb)
  Ec = dark_data.phi_b{ii}(1:Nsolid) - dark_data.phi_b{ii}(end);
  Ev = Ec-sc.Eg;
  dark_Ef{ii} = Ec(1)*ones(1,Nsolid) + Vt*log(sc.ND/sc.Nc);
  dark_Ec{ii} = Ec;
  dark_Ev{ii} = Ev;
end

% ------------------------------------------------------------------
% Compute the voltage separation between the light and dark currents
for ii=1:length(light_data.Vbb)
   Jp = light_data.Jp(ii);
   Vbb_dark(ii) = interp1(dark_data.Jp,dark_data.Vbb,Jp);
   dV_currents(ii) = Vbb_dark(ii) - light_data.Vbb(ii);
end

% ------------------------------------------------------------------
% Compute the differential capacitance, dark only
rho = {};
Q = [];
C = [];
for ii=1:length(dark_data.Vbb)
  rho = ones(1,Nsolid)*sc.ND - dark_data.n{ii} + dark_data.p{ii};
  Qsc(ii) = sum((rho(1:end-1)+rho(2:end))*.5.*diff(x_solid));
  if(ii>1)
   C(ii-1) = (Qsc(ii)-Qsc(ii-1))/(dark_data.Vbb(ii)-dark_data.Vbb(ii-1));
  end
end



% ==========================================================================================================
% Plotting

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure;   %!!!!!!!  Figure Creation  !!!!!!!
set(gcf,'Position',[1 200 500 1000])

% ------------------------------------------------------------------
% J-V Plot full range linear scale
subplot(3,1,1);
plot(light_data.Vbb,light_data.Jp/light_sys.I0,'r-',...
       dark_data.Vbb,dark_data.Jp/light_sys.I0,'b-');
ylim([min(light_data.Jp/light_sys.I0) 2])
grid on;
%legend('Photocurrent','Dark Current','Location','North');
xlabel('Bias (V)');
ylabel('Hole Current Nomalized (Arb.)');
xlim([0 2.3])
s = sprintf('Interfacial Hole Transfer Current');
%title(s);

% ------------------------------------------------------------------
% J-V plot full range log scale
subplot(3,1,2);
semilogy(light_data.Vbb,light_data.Jp/light_sys.I0,'r-',...
       dark_data.Vbb,dark_data.Jp/light_sys.I0,'b-');
ylim(0.9*[min(light_data.Jp/light_sys.I0) 4])
grid on;
%legend('Photocurrent','Dark Current','Location','South');
xlabel('Bias (V)');
ylabel('Hole Current Nomalized (Arb.)');
xlim([0 2.3])
s = sprintf('Interfacial Hole Transfer Current');
%title(s);

% ------------------------------------------------------------------
% Interfacial potential drop
subplot(3,1,3);
semilogy(light_data.Vbb,exp(0.5*light_data.Vint/Vt),'r-',...
         dark_data.Vbb,exp(0.5*dark_data.Vint/Vt),'b-');
ylim([1 100]);
grid on;
xlabel('Bias (V)');
ylabel('Bulter-Volmber Rate Contribution (1/s)');
xlim([0 2.3])
yyaxis right;
plot(light_data.Vbb,light_data.Vint,'r--',...
     dark_data.Vbb,dark_data.Vint,'b--');
ylim([0 9.1*Vt]);
ylabel('Interface Potential Drop (eV)');
s = sprintf('Butler-Volmer Kinetics');
%title(s);

%print -djpeg ResultsFig1.jpeg

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure;   %!!!!!!!  Figure Creation  !!!!!!!
set(gcf,'Position',[500 200 500 1000])

% ------------------------------------------------------------------
% J-V Plot around onset
subplot(4,1,1);
plot(light_data.Vbb,light_data.Jp/light_sys.I0,'r-');
ylim([min(light_data.Jp/max(light_data.Jp)) 1.1])
grid on;
%legend('Photocurrent','Location','South');
xlabel('Bias (V)');
ylabel('Hole Current Nomalized (Arb.)');
xlim([0 1.2])
s = sprintf('Interfacial Hole Transfer Current');
%title(s);

% ------------------------------------------------------------------
% Recombination supression plot
subplot(4,1,3)
semilogy(light_data.Vbb,light_data.R_eh_flux./light_data.G_flux,'r',...
         light_data.Vbb,light_data.G_flux./light_data.G_flux,'k')
grid on;
%legend('Recombination Flux','Generation Flux','Location','SouthWest');
xlabel('Bias (V)');
ylabel('Flux (Arb.)');
%ylim([min(light_data.R_eh_flux./light_data.G_flux) 10])
ylim([1e-5 10])
xlim([0 1.2])
s = sprintf('Depletion Generation and Recombination Flux Ratio');
%title(s);


% ------------------------------------------------------------------
% Carrier concentration plot
subplot(4,1,2)
[val,ip] = min(abs(light_data.Vbb-1.0));  % Index well beyond onset, prior to the dark current
semilogy(light_data.Vbb,light_ps/light_ps(ip),'r'); hold on;
ylim([min(light_ps/light_ps(ip))/10 2]);
%yyaxis right;
semilogy(light_data.Vbb,light_ns/max(light_ns),'k');
xlim([0 1.2])
%ylim([min(light_ps/max(light_ps)) max(light_ps/max(light_ps))*1.2]);
grid on;
ylabel('Surface Carrier Conc. (Arb.)');
xlabel('Bias (V)');
s = sprintf('Surface Carrier Concentrations');
%title(s);

% ------------------------------------------------------------------
% Quasi-Fermi level splitting plot
subplot(4,1,4)
plot(light_data.Vbb,light_dF,'k',...
     light_data.Vbb,dV_currents,'g--'); 
grid on;
xlim([0 1.2])
ylabel('\Delta V_{F} (eV)');
xlabel('Bias (V)');
ylim([0.9 1.7])
s = sprintf('Quasi Fermi Level Splitting');
%title(s);

%print -djpeg ResultsFig2.jpeg


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure; 
set(gcf,'Position',[600 100 1200 600])

% ------------------------------------------------------------------
% Illuminated band diagram plotting right after flat band
subplot(2,3,1)
[val,ip] = min(abs(light_data.Vbb-0));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,light_Ec{ip},'r'); hold on;
plot(x_solid/1e-9,light_Ev{ip},'r'); hold on;
plot(x_solid/1e-9,light_Fn{ip},'r--'); hold on;
plot(x_solid/1e-9,light_Fp{ip},'r--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(light_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Illumination Band Diagram | %1.2f V',light_data.Vbb(ip));
title(s);
hold off;


% ------------------------------------------------------------------
% Illuminated band diagram plotting just before onset
subplot(2,3,2)
[val,ip] = min(abs(light_data.Vbb-0.3));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,light_Ec{ip},'r'); hold on;
plot(x_solid/1e-9,light_Ev{ip},'r'); hold on;
plot(x_solid/1e-9,light_Fn{ip},'r--'); hold on;
plot(x_solid/1e-9,light_Fp{ip},'r--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(light_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Illumination Band Diagram | %1.2f V',light_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Illuminated band diagram plotting at onset
subplot(2,3,3)
[val,ip] = min(abs(light_data.Jp/light_sys.I0-0.5));
Vonset = light_data.Vbb(ip);
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,light_Ec{ip},'r'); hold on;
plot(x_solid/1e-9,light_Ev{ip},'r'); hold on;
plot(x_solid/1e-9,light_Fn{ip},'r--'); hold on;
plot(x_solid/1e-9,light_Fp{ip},'r--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(light_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Illumination Band Diagram | %1.2f V',light_data.Vbb(ip));
title(s);
hold off;


% ------------------------------------------------------------------
% Illuminated band diagram plotting after onset
subplot(2,3,4)
[val,ip] = min(abs(light_data.Vbb-(2*Vonset)));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,light_Ec{ip},'r'); hold on;
plot(x_solid/1e-9,light_Ev{ip},'r'); hold on;
plot(x_solid/1e-9,light_Fn{ip},'r--'); hold on;
plot(x_solid/1e-9,light_Fp{ip},'r--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(light_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Illumination Band Diagram | %1.2f V',light_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Illuminated band diagram plotting towards dark current merging
subplot(2,3,5)
[val,ip] = min(abs(light_data.Vbb-2));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,light_Ec{ip},'r'); hold on;
plot(x_solid/1e-9,light_Ev{ip},'r'); hold on;
plot(x_solid/1e-9,light_Fn{ip},'r--'); hold on;
plot(x_solid/1e-9,light_Fp{ip},'r--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(light_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Illumination Band Diagram | %1.2f V',light_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Generation term depth
G = light_sys.I0*sc.alpha*exp(-abs(x_solid)*sc.alpha);
subplot(2,3,6)
plot(x_solid/1e-9,G/(light_sys.I0*sc.alpha),'k');
ylabel('Normalized Photon Absorbtion (Arb.)');
xlabel('(nm)');

%print -djpeg ResultsFig3.jpeg


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure; 
set(gcf,'Position',[600 400 1200 600])

% ------------------------------------------------------------------
% Dark band diagram plotting right after flat band
subplot(2,3,1)
[val,ip] = min(abs(dark_data.Vbb-0));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,dark_Ec{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ev{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ef{ip},'b--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(dark_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Dark Band Diagram | %1.2f V',dark_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Dark band diagram plotting 
subplot(2,3,2)
[val,ip] = min(abs(dark_data.Vbb-0.6));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,dark_Ec{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ev{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ef{ip},'b--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(dark_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Dark Band Diagram | %1.2f V',dark_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Dark band diagram plotting 
subplot(2,3,3)
[val,ip] = min(abs(dark_data.Vbb-1.2));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,dark_Ec{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ev{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ef{ip},'b--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(dark_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Dark Band Diagram | %1.2f V',dark_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Dark band diagram plotting 
subplot(2,3,4)
[val,ip] = min(abs(dark_data.Vbb-1.8));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,dark_Ec{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ev{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ef{ip},'b--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(dark_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Dark Band Diagram | %1.2f V',dark_data.Vbb(ip));
title(s);
hold off;

% ------------------------------------------------------------------
% Dark band diagram plotting 
subplot(2,3,5)
[val,ip] = min(abs(dark_data.Vbb-2.15));
plot([-1e-15 1e-15],[1.1*min(Ev) 0],'k--'); hold on;   % Interface line
plot(x_solid/1e-9,dark_Ec{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ev{ip},'b'); hold on;
plot(x_solid/1e-9,dark_Ef{ip},'b--'); hold on;
xlim([1.01*x(1)/1e-9 x(end)/1e-9]);
ylim([1.1*min(dark_Ev{ip}) 0.1]);
xlabel('(nm)');
ylabel('With Respect to Flat Band Potential (eV)');
s = sprintf('Dark Band Diagram | %1.2f V',dark_data.Vbb(ip));
title(s);
hold off;


% ------------------------------------------------------------------
% Mott-Schottky Plot
subplot(2,3,6);
[val,ip] = min(abs(dark_data.Vbb-1.2));
invCmax2 = max(1./(C(1:ip).^2));
interpCap = linspace(0,1.0,ip);
plot(dark_data.Vbb(1:ip),(1./(C(1:ip).^2))/invCmax2,'k',...
     dark_data.Vbb(1:ip),interpCap,'k--');
ylabel('1/C^2 (Arb.)')
xlabel('Bias (V)')

%print -djpeg ResultsFig4.jpeg

%% ..................................................................
%!open R*.jpeg
