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
% o Sets default input parameters based on TaN and returns them through the system (sys), semiconductor (sc) %
%   and interface (inter) structured arrays (https://www.mathworks.com/help/matlab/ref/struct.html).         %                             
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

function [sys,sc,inter] = DefaultInputVals()

% ----------------------------------------------------------------------------------------------------------
% Fundamental Constants 
q = 1.602e-19;                           % Charge on an electron [coul] 
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m] 
kB = 1.38e-23;                           % Boltzmann constant [J/K]
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]

% ----------------------------------------------------------------------------------------------------------
% Generic system (sys) parameters, not specific to any "object"
sys.I0 = 6.2063e16*(100^2)/1;            % Incident photon flux [m^-2 s^-1]
sys.T = 300;                             % Temperature
Vt = kB*sys.T/q;                         % Thermal voltage
sys.Vbi = 0.25;                          % Built-in voltage [V], measured experimentally see: https://doi.org/10.1039/D0CP02888Fo
sys.Vmin = 0.3; %Vt; %Vt;                % Minimum bias just after flat band
sys.Vmax = 0.8;                          % Maximum bias
sys.dV = Vt/2;                           % Voltage step size

% ----------------------------------------------------------------------------------------------------------
% Calculate & set bulk semiconductor (sc) parameters
% -> Some benchmarks here... http://www.ioffe.ru/SVA/NSM/Semicond/
% -> For tantalum nitride see... https://doi.org/10.1007/s40243-016-0083-z
% -> Typical Bimolecular recombination rates are around B=1e-10 cm^3/s, see Fig. 5.4 of "Light Emitting Diodes" by F. Schubert (3rd Ed.)
mn = 2.70*m0;                                    % Bulk electron effective mass, more physically meaningful lower value when surface states are included
mp = 3.56*m0;                                    % Bulk hole effective mass, more physically meaningful lower value when surface states are included
sc.Eg = 2.2;                                   % Band gap [eV]
sc.Nc = 2*(2*pi*mn*Vt*q/(h^2))^(3/2);          % Effective density of conduction band states [#/m^3]
sc.Nv = 2*(2*pi*mp*Vt*q/(h^2))^(3/2);          % Effective density of valence band states [#/m^3]
sc.ni = sqrt(sc.Nv*sc.Nc)*exp(-sc.Eg/(2*Vt));     % Intrinsic carrier concentration [#/m^3]
sc.epsr = 17;                                     % Relative dielectric constant, see DOI: 10.1039/c4cp05616g
sc.eps = sc.epsr*eps0;                            % Overall dielectric constant
sc.alpha = 3.4e5/0.01;                            % Absorption coefficient at 400 nm wavelength [m^-1], see https://doi.org/10.1007/s40243-016-0083-z
sc.ND = 1e17*(100^3);                             % n-type: Donor density [#/m^3]
sc.p_o = (sc.ni^2)/sc.ND;                         % Bulk hole density [#/m^3]
sc.tau = 1e-8;                                    % Recombination time benchmarked to n-type bulk [s]
sc.B = 1/(sc.ND*sc.tau);                          % Bimoleular coefficient benchamarked to n-type bulk following Gartner [m^3/s] 
sc.Lp = sqrt( ((2*sc.eps)/q)*(1/sc.ND)*sc.Eg );   % Hole diffusion length [m]
sc.ratio_n = 0.99;                                % Ratio n/ND defining the limit where the electric field is considered negligible
sc.carrier_max = 1/(0.27144e-9^3);                % Density of atoms in diamond silicon [#/m^3]

% ----------------------------------------------------------------------------------------------------------
% Interface (inter) parameters
inter.kp0 = 10*1e0; %1e0 %1e1 %1e5 %1.0             % Hole transfer rate [1/s]
inter.Lint = 1e-9;                                  % Width of Butler-Volmer potential drop [m]
inter.Wp = 2.5e-9;                                  % Tunneling hole depth parameter [m] participating in surface reactions
inter.Lsc = sqrt( ((2*sc.eps)/q)*(1/sc.ND)*sc.Eg ); % Semiconductor screening length when band bending is equal to the band gap, see https://en.wikipedia.org/wiki/Depletion_region
inter.Lscreening = 10e-9; %8e-9; %30e-9; %2*5e-9;   % Hole screening length depth

% ----------------------------------------------------------------------------------------------------------
% Specialized convergence parameters for potential drop calculations
sys.beta = pi/100;                % Mixing parameter in the Gummel loop
sys.beta_phi = pi/10;             % Damping parameter for the Gummel potential guess, very sensitive convergence parameter
sys.Ngum = 200; %1000;            % Maximum number of Gummel iterations
sys.ct = 1e-4; %1e-5;             % Gummel loop convergence target
sys.Nphi = 100;                   % Maximum number of voltage guess iterations

% ========================================================================================================== %


