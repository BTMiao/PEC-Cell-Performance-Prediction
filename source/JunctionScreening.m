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
% o Returns hole and electron screening properties, as well as the semiconductor-liquid interface potential  %
%   drop (scr_data). All quantities are computed self-consistently per the light/dark conditions input.      %   
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

function [scr_data] = JunctionScreening(sys,sc,inter,grid) 

% ===========================================================================================================
% ==================....Parameters, General Variables, and Operator Definitions Start....====================

% -----------------------------------------------------------------------------------------------------------
% Fundamental parameters
k = 8.617e-5;                   % Boltzmann constant [eV/K]
eps0 = 8.85e-12;                % Permittivity of free space [F/m]
q = 1.602e-19;                  % Charge on an electron [coul]
m = 9.11e-31;                   % Electron mass (kg)
planck_constant = 6.63e-34;     % [J.s] to calculate the "effective" density of the conduction band and valence band, Nc and Nv
N_avg = 6.023e23;               % Avogadro's number (# of particle/ mole)

% -----------------------------------------------------------------------------------------------------------
% Fixed parameters, liquid properties could become input parameters in revised versions
KS_liquid = 80;               % Dielectric constant for water of the liquid  
Compensating_cation = 1.0;    % Supporting electrolyte [mol]  
Compensating_anion = 1.0;     % Supporting electrolyte [mol]  

% -----------------------------------------------------------------------------------------------------------
% General function input parameters
T = sys.T;
KS_solid = sc.epsr;
ND = sc.ND;
p_bulk = sc.p_o;
x = grid.x;
if(sys.I0 < 1e-3)
  dark_flag = 1;
else
  dark_flag = 0;
end
Vt = k*T;
kT = k*T;       % Same as Vt, different notation

% -----------------------------------------------------------------------------------------------------------
% Assign light flux parameters
G = sys.I0*sc.alpha*exp(-abs(grid.x_solid)*sc.alpha);      % Generation expression
tau = sc.tau;
Lint =  inter.Lint;
Wp = inter.Wp;
kp0 = inter.kp0;
p_eq = grid.p_eq;
if(max(abs(p_eq))==0)    % Pertaining to Built-in voltage calculations
  dark_flag = 1;
end
[val,BV_ii_start] = min(abs(x+Lint/2));    % Defining the Butler-Volmer potential drop contribution to kinetics
[val,BV_ii_end] = min(abs(x-Lint/2));      % Defining the Butler-Volmer potential drop contribution to kinetics

% -----------------------------------------------------------------------------------------------------------
% Spatial grid defintion for the semiconductor and liquid, depends on system screening, all lengths in meters.
[val,Nsolid] = min(abs(x));             % Number of grid points inside solid
N = length(x);                          % Total number of grid point  
Nliquid = N-Nsolid;                     % Number of grid points inside liquid
x_solid = x(1:Nsolid);                   % Semiconductor region
h = zeros(1,N+1);                         
hd = diff(x);                           % Grid spacing
h = [hd(1) hd hd(end)];                 % Utilized in the Laplacian operator

% -----------------------------------------------------------------------------------------------------------
% Variable dielectric constant distribution
KSmat = x*0;
KSmat(1:Nsolid) = KS_solid;            % Dielectric constant matrix of semiconductor 
KSmat(Nsolid+1:N) = KS_liquid;         % Dielectric constant matrix of liquid
deps = 2.5e-11;                        % Making sure that the dielectric constant changes smoothly at the interface to avoid numerical error. 2.5e-11 m results in a 2 Ang transtion region between dielectric constants 
feps = 1./(1+exp(-x/deps));            % Assumes that the solid dielectric constant is less than the liquid
KSmat = feps*(KS_liquid - KS_solid) + KS_solid;   
eps = KSmat*eps0;                      % Overall variable dielectric constant in F/m

% -----------------------------------------------------------------------------------------------------------
% Modified Laplacian for variable dielectric
La = (2./(h(2:N).*(h(2:N)+h(3:N+1)))); % Lower Laplacian diagonal                                         
Lb = (2./(h(1:N).*h(2:N+1)));          % Centre Laplacian diagonal
Lc = 2./(h(2:N).*(h(1:N-1)+h(2:N)));   % Upper Laplacian diagonal
laplacian_dirichlet1 = -diag(Lb) + diag(Lc,1) + diag(La,-1);
laplacian_dirichlet = diag(eps)*laplacian_dirichlet1;

% -----------------------------------------------------------------------------------------------------------
% Gradient Operator
Ga = [0 1./(h(3:N)) 0];
Gc = [0 1./(h(3:N))];
grad_op = -diag(Ga) + diag(Gc,1);       
grad_coeff = (grad_op*eps')';
Ga_eps = [0 grad_coeff(2:N-1).*(1./(h(3:N))) 0];
Gc_eps = [0 grad_coeff(2:N-1).*(1./(h(3:N)))];
grad_eps = -diag(Ga_eps) + diag(Gc_eps,1);

% -----------------------------------------------------------------------------------------------------------
% Semiconductor charge parameters
accep = x*0;                     % n-type semiconductor assumed
donor = x*0;
donor(1:Nsolid) = ND;            %  Donor concentration in the semiconductor region
nuclear_static = accep + donor;  % All unscreened nuclear charges

% ====================....Parameters, General Variables, and Operator Definitions End....====================
% ===========================================================================================================


% ===========================================================================================================
% =========================....Potential Profile & Convergence Variables Start....===========================

% -----------------------------------------------------------------------------------------------------------
% Convergence parameters
Ngum = sys.Ngum;
beta = sys.beta;       % Gummel loop mixing parameter
beta_start = beta;
beta_phi = sys.beta_phi;   % Convergence can be very sensitive to this parameter
beta_phi_start = beta_phi;
ct = sys.ct; %1e-4; %1e-3; %kT/2500; %5e-3; %kT/2500;

% -----------------------------------------------------------------------------------------------------------
% Spatially dependent screening quantities, for charge conservation checks
electrode_charge = zeros(1,N); % Electrode charge density at a given bias (includes ss charge)
liquid_charge = zeros(1,N);    % Liquid charge density at a given bias
electrode_total = zeros(1,1);  % Electrode charge at a given bias (includes ss charge)
liquid_total = zeros(1,1);     % Liquid charge at a given bias
                               % positive bias lowers the semiconductor Fermi level down 
                               % relative to the liquid Fermi level (and vice versa)
   
% ==========================....Potential Profile & Convergence Variables End....============================
% ===========================================================================================================


% ===========================================================================================================
% ==================....Potential Profile (Zero and Applied Bias) Calculations Start....=====================
% Start with intial potential guess
phi_guess = grid.phi_guess;
dphi_old = zeros(N,1);
p_max = sc.carrier_max;    % Hole density upper bound, e.g. based on density of atoms in Si diamond lattice
fileID = fopen(sys.RunDataFile,'a');  % Storing the convergence data through the append ('a') option

% -----------------------------------------------------------------------------------------------------------
% Begin the bias calculation
tol = 1e10;                            % Inner Gummel loop SCF tolerance result
tolViter = 1e10;                       % Outer potential convergence (Gummel) tolerance result
phi = -1*(phi_guess-phi_guess(end));   % Initial guess 
phi_old_Viter = phi;                   % Old initial guess for the outer potential loop
Vbias = phi_guess(end)-phi_guess(1);   % Applied bias
gg = 1;                                % Number of Gummel loops

while((tolViter>ct) && (gg<sys.Nphi))  % Outer potential SCF loop
   
   % --------------------------------------------------------------------------------------------------------
   % Begin the Gummel loop
   ii = 1;
   while((tol>ct) && (ii<sys.Ngum))  % Gummel loop

      % .....................................................................................................
      % Semiconductor depletion region setup
      n = [ND*exp(+(phi(1:Nsolid)-phi(1))/kT) zeros(1,Nliquid)];   % Electron concentration
      if(n(Nsolid)>=sc.ratio_n*ND)
        xd_ii = Nsolid-1;
        xd_iiminus1 = xd_ii-1;
        xd_pos = (x_solid(xd_ii) + x_solid(xd_ii-1))/2;  % Place in the middle
      else  % Positioning of the depletion region width
        [val,xd_ii] = min(abs(n(1:Nsolid)-sc.ratio_n*ND));
        if(xd_ii==1)
          xd_ii=2;
        end
        temp_n = n(1:Nsolid)+linspace(1/Nsolid,1,Nsolid)*ND*1e-10;   % Provides a unique n array for interpolation
        xd_pos = interp1(temp_n(xd_ii-1:Nsolid),x_solid(xd_ii-1:Nsolid),sc.ratio_n*ND);
        if(x_solid(xd_ii) < xd_pos)  % Find the relative position between two grid points
          xd_ii = xd_ii + 1;
        end
        xd_iiminus1 =  xd_ii-1;
      end
      zeta = (x_solid(xd_ii)-xd_pos)/(x_solid(xd_ii)-x_solid(xd_iiminus1));

      % .....................................................................................................
      % Semiconductor screening computation in the depletion region
      if(dark_flag > 0) % Dark condition quantities
        p = [p_bulk*exp(-(phi(1:Nsolid)-phi(1))/kT) zeros(1,Nliquid)];   % hole concentration
        p_s = p(Nsolid);
      else  % Compute the hole density using the illuminated quantities
        G_fluxA = sum((G(xd_ii:end-1)+G(xd_ii+1:end))*.5.*diff(x_solid(xd_ii:Nsolid)));
        G_fluxB = sum((G(xd_iiminus1:end-1)+G(xd_iiminus1+1:end))*.5.*diff(x_solid(xd_iiminus1:Nsolid)));
        G_flux = (1-zeta)*G_fluxA+zeta*G_fluxB;
        Ec = -phi_old_Viter(1:Nsolid);
        Ec = Ec - Ec(1);   
        Vint = abs(phi_old_Viter(BV_ii_end)-phi_old_Viter(BV_ii_start));  % Interface potential term
        numerator_term1 = G_flux + p_bulk*abs(xd_pos)/tau;
        num_term2 = p_eq.*exp(-abs(x_solid)/Wp)*kp0*exp(Vint/(2*Vt));
        numerator_term2A = sum((num_term2(xd_ii:end-1)+num_term2(xd_ii+1:end))*.5.*diff(x_solid(xd_ii:end)));
        numerator_term2B = sum((num_term2(xd_iiminus1:end-1)+num_term2(xd_iiminus1+1:end))*.5.*diff(x_solid(xd_iiminus1:Nsolid)));
        numerator_term2 = (1-zeta)*numerator_term2A+zeta*numerator_term2B;
        numerator = numerator_term1 + numerator_term2;
        denom_term1 = (n(1:Nsolid).*exp((Ec-Ec(end))/Vt))/(ND*tau);
        denom_term2 = kp0*exp(Vint/(2*Vt))*(exp((Ec-Ec(end))/Vt).*exp(-abs(x_solid)/Wp));
        denominatorA = sum((denom_term1(xd_ii:end-1)+denom_term1(xd_ii+1:end))*.5.*diff(x_solid(xd_ii:end))) ...
                   + sum((denom_term2(xd_ii:end-1)+denom_term2(xd_ii+1:end))*.5.*diff(x_solid(xd_ii:end)));
        denominatorB = sum((denom_term1(xd_iiminus1:end-1)+denom_term1(xd_iiminus1+1:end))*.5.*diff(x_solid(xd_iiminus1:end))) ...
                   + sum((denom_term2(xd_iiminus1:end-1)+denom_term2(xd_iiminus1+1:end))*.5.*diff(x_solid(xd_iiminus1:end)));
        denominator = (1-zeta)*denominatorA+zeta*denominatorB;
        p_s = numerator/denominator;
        p_dark = [p_bulk*exp(-(phi(1:Nsolid)-phi(1))/kT) zeros(1,Nliquid)];   % hole concentration
        if(p_dark(Nsolid) > 0.001*p_s)   % Gradually phase in the dark quasi-Fermi level
          p_s = p_s+p_dark(Nsolid);
        end
        p = 0*x;
        p(1:Nsolid) = p_s*exp(-(phi(1:Nsolid)-phi_old_Viter(Nsolid))/Vt);   % Need to use stable reference potential for convergence
      end

      % .....................................................................................................
      % Semiconductor hole screening properties beyond the depletion region and overall screening limits
      if(p(Nsolid) > p_max)   % Sets a "physical" upper bound on the allowed hole density, also needed for quasi-Fermi level convergence at high biases
        p_s = p_max;
        p(1:Nsolid) = p_s*exp(-(phi(1:Nsolid)-phi(Nsolid))/Vt);
      end
      if(~(dark_flag>0))   % Compute the hole concentration beyond the depletion region
        scr_data.p = p;   % Intial flat everywhere hole quasi Fermi level guess
        scr_data.n = n;
        scr_data.p_s = p_s;
        xd_ii_temp = xd_ii;
        scr_data.xd_ii = xd_ii_temp;
        p_fb = HoleFlatBandConc(sys,sc,inter,grid,scr_data);
        Nfb = length(p_fb);
        p(1:Nfb) = p_fb;
      end 
      nuclear = nuclear_static;
      solid_charge=p-n+nuclear;
      %%%solid_charge=0*p-n+nuclear;   % Turns off hole inversion screening

      % .....................................................................................................
      % Solution concentration perturbations from the bulk (requires a high supporting electrolyte concentration)!!
      phi_liq_bulk = phi(end);
      cation = [zeros(1,Nsolid) ...     % Supportiong cation concentration
                Compensating_cation*1000*N_avg*exp(-(phi(Nsolid+1:N)-phi_liq_bulk)/kT)];     
      anion = [zeros(1,Nsolid) ...      % Supportiong anion concentration
               Compensating_anion*1000*N_avg*exp((phi(Nsolid+1:N)-phi_liq_bulk)/kT)];

      % .....................................................................................................
      % LU decomposition solution to the potential
      A = laplacian_dirichlet+1*grad_eps - diag(q*(cation+anion+p+n)/kT);                                                                  
      %%%A = laplacian_dirichlet+1*grad_eps - diag(q*(cation+anion+0*p+n)/kT);  % Turns off hole inversion screening                                                                  
      A(1,:) = [1 zeros(1,N-1)];   % Use the Dirichlet boundary condition (semiconductor bulk) 
      A(N,:) = [zeros(1,N-1) 1];   % Use the Dirichlet boundary condition (liquid bulk)
      vector = (-1*laplacian_dirichlet*(phi'))-1*grad_eps*(phi') + (q*(anion-cation-solid_charge))';
      vector(1) = 0;               % Dirichlet boundary condition continued
      vector(end) = 0;             % Dirichlet boundary condition continued
      dphi = (LUdec(diag(A,-1),diag(A),diag(A,1),vector))';% LU decomposition 
      phi = phi+beta*(dphi');                                       
      tol = max(abs(dphi));
      printstring = sprintf('-> Bias: %f | gg: %i | ii: %i | tol: %f | tolViter: %f | beta: %f | beta_phi: %f',....
                            Vbias,gg,ii,tol,tolViter,beta,beta_phi);
      disp(printstring);
      printstring = cat(2,printstring,'\n');
      fprintf(fileID,printstring);

      if((ii<10) && (tol<ct))  % Execute at least 10 Gummel iterations
        tol = 10*ct;
      end
      ii = ii + 1;
   end % Gummel loop termination


   % -----------------------------------------------------------------------------------
   % Reduce the mixing parameter if the maxium number of Gummel iterations is reached
   if(ii>(Ngum-1))
     if(beta > 4*ct)
       beta = beta/2;
     end
   else
     if(beta<(beta_start/2))
       beta = beta*2;
     end
   end

   % -----------------------------------------------------------------------------------
   % Refine the reference potential guess
   tolViter = max(abs(phi_old_Viter-phi));   
   phi_old_Viter = (1-beta_phi)*phi_old_Viter + beta_phi*phi;
   if( (gg > 2) && (ii>(Ngum-1)) )
     if(beta_phi > 4*ct)
       beta_phi = beta_phi/2;
     end
   elseif( (gg > 2) && (beta_phi<(beta_phi_start/2)) )
     beta_phi = beta_phi*2;
   end
   gg = gg + 1;
   
   printstring = sprintf('-> Bias: %f | gg: %i | ii: %i | tol: %f | tolViter: %f | beta: %f | beta_phi: %f',....
                            Vbias,gg,ii,tol,tolViter,beta,beta_phi);
   disp(printstring);
   printstring = cat(2,printstring,'\n');
   fprintf(fileID,printstring);

end % Overall potential profile convergence loop

% --------------------------------------------------------------------------------------------------------
% Post processing calculation results & checks
phi = -1*(phi-phi(1));               % Band bending with the semiconductor bulk at 0 eV
electrode_charge = p-n+nuclear;      % Charge on the electrode
liquid_charge = cation-anion;        % Charge on the electrode    
electrode_total = q*sum((electrode_charge(1:end-1)+electrode_charge(2:end))*.5.*diff(x));
liquid_total = q*sum((liquid_charge(1:end-1)+liquid_charge(2:end))*.5.*diff(x));

% --------------------------------------------------------------------------------------------------------
% Assign data to return and file closing
scr_data.phi = phi;
scr_data.p_s = p_s;
scr_data.n = n;
scr_data.p = p;
scr_data.xd_ii = xd_ii;
scr_data.xd_pos = xd_pos;
fclose(fileID);

% ===================....Potential Profile (Zero and Applied Bias) Calculations End....======================
% ===========================================================================================================

