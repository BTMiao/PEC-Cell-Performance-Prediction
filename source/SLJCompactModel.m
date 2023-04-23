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
% o Primary execution function for the model, takes only the input parameter file and calls all other        %
%   functions to run the simulation and then save the final output data.                                     %
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

function [] = SLJCompactModel(inputfile);    % Saves data to file rather than returning it

% ----------------------------------------------------------------------------------------------------------
% Fundamental Constants 
q = 1.602e-19;                           % Charge on an electron [coul] 
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m] 
kB = 1.38e-23;                           % Boltzmann constant [J/K]
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]

% ----------------------------------------------------------------------------------------------------------
[sys,sc,inter] = ParseInputFile(inputfile);

% ----------------------------------------------------------------------------------------------------------
% Setup the grid and space charge region, needed for plotting
xgrid = GetSpatialGrid(sys,sc,inter);

% ----------------------------------------------------------------------------------------------------------
% Get the output convergence run data ready for output
fileID = fopen(sys.RunDataFile,'w');
fprintf(fileID,'');
fclose(fileID);

% ----------------------------------------------------------------------------------------------------------
% Execute the run
tic;
data = SLJSemiAnalytical(sys,sc,inter,xgrid);
data.run_time = toc;
fileID = fopen(sys.RunDataFile,'a');
fprintf(fileID,'Calculation Run Time (s) %f | Light Flux: %f\n\n',data.run_time,sys.I0);
fclose(fileID);

% ----------------------------------------------------------------------------------------------------------
% Formatting output variables  (more could be returned, such as quasi-Fermi levels)
RunDataFile = sys.RunDataFile;
jv_data.Jp = q*data.Jp;
jv_data.G_flux = data.G_flux;
jv_data.R_eh_flux = data.R_eh_flux;
jv_data.R_eh = data.R_eh;
jv_data.n = data.n;
jv_data.p = data.p;
jv_data.phi_b = data.phi_b;
jv_data.Vbb = data.Vbb;
jv_data.Vint = data.Vint;
jv_data.phi_b_eq = data.phi_b_eq;
jv_data.p_eq = data.p_eq;
jv_data.run_time = data.run_time;
jv_data.x = xgrid.x;
inter_temp.kp0 = inter.kp0;
inter_temp.Lint = inter.Lint;
inter_temp.Wp = inter.Wp;
inter_temp.Lscreening = inter.Lscreening;
inter = inter_temp;
sc_temp.Eg = sc.Eg;
sc_temp.Nc = sc.Nc;
sc_temp.Nv = sc.Nv;
sc_temp.epsr = sc.epsr;
sc_temp.alpha = sc.alpha;
sc_temp.ND = sc.ND;
sc_temp.tau = sc.tau;
sc_temp.ratio_n = sc.ratio_n;
sc_temp.carrier_max = sc.carrier_max;
sc = sc_temp;
sys_temp.I0 = sys.I0;
sys_temp.T = sys.T;
sys_temp.Vmin = sys.Vmin; 
sys_temp.Vmax = sys.Vmax;
sys_temp.Vbi = sys.Vbi;
sys_temp.dV = sys.dV;
sys_temp.beta = sys.beta;
sys_temp.beta_phi = sys.beta_phi;
sys_temp.Ngum = sys.Ngum;
sys_temp.ct = sys.ct; 
sys_temp.Nphi = sys.Nphi;
sys = sys_temp;

% ----------------------------------------------------------------------------------------------------------
% File writing
names = split(RunDataFile,'.');
filename = cat(2,names{1},'.mat');
save(filename,'jv_data','sc','sys','inter','q','kB','eps0','h','m0');   % Saves all results to a single Matlab file for plotting and processing

% ========================================================================================================== %


