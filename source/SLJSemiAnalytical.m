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
% o Computes and returns all simulation data pertaining to the bias range input and the illumination         %
%   intensity provided.  The voltage sweep loop is implemented in this function.                             %
% o Takes as input the system (sys), semiconductor (sc), interface (inter) parameters as defined by the      %
%   input file via ParseInputFile.m through or directly via DefaultInputVals.m.  It also takes input the     %
%   spatial grid (grid) assigned in GetSpatialGrid.m.                                                        %
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

function [data] = SLJSemiAnalytical(sys,sc,inter,grid) 

% ----------------------------------------------------------------------------------------------------------
% Fundamental Constants & Parameters
q = 1.602e-19;                           % Charge on an electron [coul] 
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m] 
kB = 1.38e-23;                           % Boltzmann constant [J/K]
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]
Vt = kB*sys.T/q;                         % Thermal voltage

% ----------------------------------------------------------------------------------------------------------
% Assign shorthand internal function grid variables
Nsolid = grid.Nsolid;
N = grid.N;
hd = grid.hd;

% ----------------------------------------------------------------------------------------------------------
% Illumination quantities
G = sys.I0*sc.alpha*exp(-abs(grid.x_solid)*sc.alpha);      % Generation expression
G_flux_tot = sum((G(1:end-1)+G(2:end))*.5.*diff(grid.x_solid));

% ----------------------------------------------------------------------------------------------------------
% Compute equilibrium screening quantities, using an analytical guess for the potential
fileID = fopen(sys.RunDataFile,'a');  % Storing the convergence data through the append ('a') option
fprintf(fileID,'************* START Built-in Voltage Convergence Results ******************\n');
fclose(fileID);
xd_eq = sqrt((2*sc.eps*(sys.Vbi))/(q*sc.ND));    
phi_guess = zeros(1,grid.N);
phi_guess = ((q*sc.ND)/(2*sc.eps))*((grid.x+xd_eq).^2);  
[val,xd_index] = min(abs(grid.x+xd_eq));
if(grid.x(xd_index)<=-xd_eq)
  phi_guess(1:xd_index) = 0;
else
  phi_guess(1:xd_index-1) = 0;
end  
phi_guess(Nsolid+1:end) = phi_guess(Nsolid);
phi_guess = phi_guess - phi_guess(1);                       % Zero at the contact electrode
grid.phi_guess = (sys.Vbi/(phi_guess(end)-phi_guess(1)))*phi_guess;  % Scale to correct Vbi (redundant check)
grid.p_eq = zeros(1,Nsolid);  % Only zero for dark calculations
scr_data = JunctionScreening(sys,sc,inter,grid);
phi_eq = scr_data.phi;
Ec_eq = phi_eq(1:Nsolid);              %  Center 0eV at the top of the conduction band minimum (bulk)
Ev_eq = phi_eq(1:Nsolid)-sc.Eg;           %  Valance band shifted 
Ef_eq = Ec_eq(1) + Vt*log(sc.ND/sc.Nc);
n_eq = sc.ND*exp(-phi_eq(1:Nsolid)/Vt);
p_eq = (sc.ni^2)./n_eq;
grid.p_eq = p_eq;    % Needed for all later illumination calculations
fileID = fopen(sys.RunDataFile,'a');  % Storing the convergence data through the append ('a') option
fprintf(fileID,'************* END Built-in Voltage Convergence Results ******************\n');
fclose(fileID);

% ----------------------------------------------------------------------------------------------------------
% Compute Bias Quantities (either dark or under illumination)
fileID = fopen(sys.RunDataFile,'a');  % Storing the convergence data through the append ('a') option
if(sys.I0 < 1e-3)
  fprintf(fileID,'************* START Dark Current Convergence Results ******************\n');
else
  fprintf(fileID,'************* START Photocurrent Convergence Results ******************\n');
end
fclose(fileID);
Vbb_min = sys.Vmin;
Vbb_max = sys.Vmax; 
dVbb = sys.dV; 
Vbb = Vbb_min;
Vbb_old = Vbb_min-dVbb;   % Sets first voltage at Vbb_min
phi = zeros(1,N);           % Total potential drop
iv = 1;
while(Vbb<=Vbb_max)   % Start of the bias loop 

  % ------------------------------------------------------------------------------------------
  % Bias step clamping for the onset region, the clamping parameters are quite sensitive wrt dVbb
  if(Vbb<=(Vbb_min+dVbb))  % Use an analytical guess  
    phi_guess = phi_eq; 
    phi_guess = (Vbb/(phi_guess(end)-phi_guess(1)))*phi_guess;
    grid.phi_guess = phi_guess;    
  else   % Use the old result as a potential guess
    phi_guess = phi*(Vbb/Vbb_old);
    phi_guess = phi_guess - phi_guess(1);
    grid.phi_guess = phi_guess;    
  end 
  
  [scr_data] = JunctionScreening(sys,sc,inter,grid);
  grid.phi_guess = scr_data.phi;
  [scr_data] = JunctionScreening(sys,sc,inter,grid);
  phi = scr_data.phi;
  Vint = abs(phi(grid.BV_ii_end)-phi(grid.BV_ii_start));
  Ec = phi(1:Nsolid);              %  Center 0eV at the top of the conduction band minimum (bulk)
  Ev = phi(1:Nsolid)-sc.Eg;        %  Valance band shifted
  n = scr_data.n(1:Nsolid);
  xd_ii = scr_data.xd_ii;
  xd_iiminus1 = xd_ii-1;
  xd_pos = scr_data.xd_pos;
  zeta = (grid.x_solid(xd_ii)-xd_pos)/(grid.x_solid(xd_ii)-grid.x_solid(xd_iiminus1));
  p_s = scr_data.p_s;
  p = scr_data.p(1:Nsolid);
  Fn = Ec(1) + Vt*log(sc.ND/sc.Nc);
  dFp = -Vt*log(p_s/sc.Nv);
  Fp_shift = (sc.Eg-abs(Fn)) - [abs(phi(Nsolid)-phi(1))]-dFp;
  Fp = Fn-Fp_shift*ones(1,Nsolid);

  if(sys.I0 < 1e-3)  % Dark current calculation
    R_eh = zeros(1,Nsolid);
    R_ht = inter.kp0*exp(+0.5*Vint/Vt)*[(p-p_eq).*exp(-abs(grid.x_solid)/inter.Wp)];     % Hole transfer term
    G_flux = 0;
  else  % Iluminated junction calculation
    if(xd_ii == Nsolid) 
      Nfb = xd_ii;    % Number of points in the flat-band (fb) region
    else
      Nfb = xd_ii-1;  % Enforce zero diffusion current between two regions
    end
    Fp(1:Nfb) = Ev(1:Nfb)-Vt*log(p(1:Nfb)/sc.Nv);
    R_eh = ((n.*p)/sc.ND-sc.p_o)/sc.tau;
    G_fluxA = sum((G(xd_ii:end-1)+G(xd_ii+1:end))*.5.*diff(grid.x_solid(xd_ii:Nsolid)));
    G_fluxB = sum((G(xd_iiminus1:end-1)+G(xd_iiminus1+1:end))*.5.*diff(grid.x_solid(xd_iiminus1:Nsolid)));
    G_flux = (1-zeta)*G_fluxA+zeta*G_fluxB;
  end   % End illumination calculation
  R_ht = inter.kp0*exp(+0.5*Vint/Vt)*[(p-p_eq).*exp(-abs(grid.x_solid)/inter.Wp)];     % Hole transfer term
  R = R_eh + R_ht;

  % ------------------------------------------------------------------------------------------
  % Store all quantities
  Jp_list(iv) = sum((R_ht(1:end-1)+R_ht(2:end))*.5.*diff(grid.x_solid));   % Total hole current
  G_flux_list(iv) = G_flux;  % Varies when flux is contained to the depletion region
  R_eh_fluxA = sum((R_eh(xd_ii:end-1)+R_eh(xd_ii+1:end))*.5.*diff(grid.x_solid(xd_ii:Nsolid)));
  R_eh_fluxB = sum((R_eh(xd_iiminus1:end-1)+R_eh(xd_iiminus1+1:end))*.5.*diff(grid.x_solid(xd_iiminus1:Nsolid)));
  R_eh_flux_list(iv) = (1-zeta)*R_eh_fluxA+zeta*R_eh_fluxB;
  R_eh_list{iv} = R_eh;
  dF_list(iv) = Fp_shift;  
  p_s_list(iv) = p_s;
  n_list{iv} = n;
  p_list{iv} = p;
  Fn_list(iv) = Fn;
  Fp_list{iv} = Fp;
  Ev_list{iv} = Ev;
  Ec_list{iv} = Ec;
  Bias(iv) = abs(phi(end)-phi(1));
  phi_list{iv} = phi;
  dphi_list{iv} = phi-phi_guess;
  Vint_list(iv) = Vint;

  iv = iv+1;
  Vbb = Vbb_min + iv*dVbb;
  Vbb_old = Bias(iv-1);
end   % Conclusion of the bias loop
NV = iv-1;   % Number of voltage steps


% ----------------------------------------------------------------------------------------------------------
% Return the data
data.Jp = Jp_list;
data.G_flux = G_flux_list;  % Varies when flux is contained to the depletion region
data.R_eh_flux = R_eh_flux_list;
data.R_eh = R_eh_list;
data.dF = dF_list;
data.p_s = p_s_list;
data.n = n_list;
data.p = p_list;
data.Fn = Fn_list;
data.Fp = Fp_list;
data.Ev = Ev_list;
data.Ec = Ec_list;
data.Vbb = Bias;
data.phi_b = phi_list;
data.Vint = Vint_list;
data.NV = NV;
data.phi_b_eq = phi_eq;  % Built-in voltage potential profile
data.p_eq = p_eq;

% ----------------------------------------------------------------------------------------------------------
% Write data
fileID = fopen(sys.RunDataFile,'a');  % Storing the convergence data through the append ('a') option
if(sys.I0 < 1e-3)
  fprintf(fileID,'************* END Dark Current Convergence Results ******************\n');
else
  fprintf(fileID,'************* END Photocurrent Convergence Results ******************\n');
end
fclose(fileID);

% ========================================================================================================== %

