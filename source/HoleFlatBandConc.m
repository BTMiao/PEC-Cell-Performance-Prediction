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
% o Returns the hole concentration beyound the depletion region where the semiconductor bands are flat.      %
% o Takes as input the system (sys), semiconductor (sc), interface (inter) parameters as defined by the      %
%   input file via ParseInputFile.m through or directly via DefaultInputVals.m.  It also takes input the     %
%   spatial grid (grid) assigned in GetSpatialGrid.m, and screening data (scr_data) determined by            %
%   JunctionScreening.m.                                                                                     %
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

function [p_fb] = HoleFlatBandConc(sys,sc,inter,grid,scr_data) 

% ----------------------------------------------------------------------------------------------------------
% Fundamental Constants & Parameters
q = 1.602e-19;                           % Charge on an electron [coul] 
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m] 
kB = 1.38e-23;                           % Boltzmann constant [J/K]
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]
Vt = kB*sys.T/q;                         % Thermal voltage

% ----------------------------------------------------------------------------------------------------------
% Function shorthand input paramters, used internally
Nsolid = length(grid.x_solid);
xd_ii = scr_data.xd_ii;
Lp2 = (sc.Lp)^2;
n = scr_data.n;
p = scr_data.p;
G = sys.I0*sc.alpha*exp(-abs(grid.x_solid)*sc.alpha);  % Usually neglect the generation term, as at high bias it can provide a non-intuitive perturbed hole quasi-Fermi level in the flat-band region
hd = grid.hd;
p_eq = grid.p_eq;
dx = hd(end);

% ----------------------------------------------------------------------------------------------------------
% Only use the flux accounting for 99% of the incoming light to avoid quasi-Fermi level descrepancies at
% very large forward biases...
G_flux_tot = sum((G(1:end-1)+G(2:end))*.5.*diff(grid.x_solid(1:Nsolid)));
G_flux_frac = sum((G(xd_ii:end-1)+G(xd_ii+1:end))*.5.*diff(grid.x_solid(xd_ii:Nsolid)));
ratio = G_flux_frac/G_flux_tot;
jj = 0;
while(ratio < 0.99)   % The remaining 1% of flux is neglible with respect to the overall current
   jj = jj+1;
   G_flux_frac = sum((G(xd_ii-jj:end-1)+G(xd_ii+1-jj:end))*.5.*diff(grid.x_solid(xd_ii-jj:Nsolid)));
   ratio = G_flux_frac/G_flux_tot;
end
G(1:xd_ii-jj) = 0; % Ignore the flux tail

% ----------------------------------------------------------------------------------------------------------
% START Diffusion operator calculation
if(xd_ii == Nsolid) 
  Nfb = xd_ii;    % Number of points in the flat-band (fb) region
else
  Nfb = xd_ii-1;  % Enforce zero diffusion current between two regions
end

Za = -Lp2*(2./(hd(2:Nfb).*(hd(2:Nfb)+hd(3:Nfb+1)))); % Lower Laplacian diagonal                                         
Zb = Lp2*(2./(hd(1:Nfb).*hd(2:Nfb+1)));          % Centre Laplacian diagonal
Zc = -Lp2*(2./(hd(2:Nfb).*(hd(1:Nfb-1)+hd(2:Nfb))));   % Upper Laplacian diagonal
 
% ..........................................................................................
% Hole density dependent terms
%Zb = Zb + 1*ones(1,Nfb) + sc.tau*inter.kp0*exp(-abs(x(1:Nfb)/Wp);  % Diagonal, Gartner recombination term 
Zb = Zb + (n(1:Nfb)/sc.ND) + sc.tau*inter.kp0*exp(-abs(grid.x(1:Nfb))/inter.Wp);  % Diagonal, Binomial recombination term 
  
% ..........................................................................................
% Compute the hole distribution in the flat-band region
Up = sc.tau*(G(1:Nfb)') + ones(Nfb,1)*sc.p_o ...
      + sc.tau*inter.kp0*p_eq(1:Nfb)*(exp(-abs(grid.x(1:Nfb))/inter.Wp)');  % Should be a single column vector
if(xd_ii == Nsolid)   % Neumann BC
  Zb(1) = 1;
  Zc(1) = 0;
  Zb(end) = Zb(end)-1*((Lp2)/(dx^2));   % Will need to change when variable grid spacing is utilized
  Up(1) = sc.p_o;
  %disp('Neumann BC called');
else
  Zb(1) = 1;
  Zc(1) = 0;
  Zb(end) = 1;
  Za(end) = 0;
  Up(1) = sc.p_o;
  Up(end) = p(xd_ii);  % set by depletion region quasi-Fermi level position
  %disp('Dirichlet BC called');
end
p_fb = LUdec(Za,Zb,Zc,Up);

% END Diffusion operator calculation
% ----------------------------------------------------------------------------------------------------------


% ========================================================================================================== %


