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
% o Returns the spatial grid for simulating the semiconductor-liquid junction.                               %
% o Takes as input the system (sys), semiconductor (sc), interface (inter) parameters as defined by the      %
%   input file via ParseInputFile.m through or directly via DefaultInputVals.m.                              %       
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

function [grid] = GetSpatialGrid(sys,sc,inter) 

% -----------------------------------------------------------------------------------------------------------
% Fundamental parameters
k = 8.617e-5;                   % Boltzmann constant [eV/K]
eps0 = 8.85e-12;                % Permittivity of free space [F/m]
q = 1.602e-19;                  % Charge on an electron [coul]
m = 9.11e-31;                   % Electron mass (kg)
planck_constant = 6.63e-34;     % [J.s] to calculate the "effective" density of the conduction band and valence band, Nc and Nv
N_avg = 6.023e23;               % Avogadro's number (# of particle/ mole)

% -----------------------------------------------------------------------------------------------------------
% Function shorthand input paramters, used internally
Eg = sc.Eg;
eps = sc.eps;
ND = sc.ND;
Lscreening = inter.Lscreening;   % Assumed hole inversion screening length
Lint = inter.Lint;               % Assumed hole inversion screening length

% -----------------------------------------------------------------------------------------------------------
% Setup the grid and space charge region
Lsc = sqrt( ((2*eps)/q)*(1/ND)*Eg );   % Semiconductor screening length when band bending is equal to the band gap, see https://en.wikipedia.org/wiki/Depletion_region
dx_1 = Lsc/100; %/100                  % This parameter should not be too coarse
Nsr = round((Lscreening/0.25e-10)/100)*100+1;
x_1 = [-3*Lsc-Lscreening:dx_1:-Lscreening];
vals = [x_1(1)-dx_1:dx_1/10:x_1(1)];
x_1 = cat(2,vals(1:end-1),x_1);               % Hole quasi-Fermi level drop region
x_2 = linspace(-Lscreening,0,Nsr);            % Define grid points in the high electric field region of the semiconductor
x_solid = [x_1(1:end) x_2(2:end)];            % Grid points in the "solid" region
dx_2 = x_2(2)-x_2(1);
x_3 = [0:dx_2:Lscreening];                    % Define grid points in the Helmholtz layer of the liquid
x_4 = linspace(x_3(end),x_3(end)+0.1e-7,81);  % Define grid points in the bulk of the liquid 
x_liquid = [x_3 x_4(2:end)];                  % Grid points at the liquid side
x = [x_solid x_liquid(2:end)];                % Grid points in total

% -----------------------------------------------------------------------------------------------------------
% Butler-Volmer and grid difference parameters
[val,BV_ii_start] = min(abs(x+Lint/2));    % Defining the Butler-Volmer potential drop contribution to kinetics
[val,BV_ii_end] = min(abs(x-Lint/2));      % Defining the Butler-Volmer potential drop contribution to kinetics
N = length(x);                             % Total number of grid points
[val,Nsolid] = min(abs(x));                % Number of grid points in the solid
x_solid = x(1:Nsolid);                     % Grid in the semiconductor/solid region
hd = diff(x_solid);                        % Grid spacing in the continuity equation Laplacian
hd = [hd(1) hd hd(end)];

% -----------------------------------------------------------------------------------------------------------
% Assign "object" quantities for returning
grid.x = x;
grid.N = N;
grid.x_solid = x_solid;
grid.Nsolid = Nsolid;
grid.BV_ii_start = BV_ii_start;
grid.BV_ii_end = BV_ii_end;
grid.hd = hd;

% ========================================================================================================== %



