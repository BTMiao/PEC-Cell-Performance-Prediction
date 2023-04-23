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
% o Parses an input file into the structured arrays: system (sys), semiconductor (sc) & interface (inter).   %
%   https://www.mathworks.com/help/matlab/ref/struct.html                                                    %
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

function [sys,sc,inter] = ParleInputFile(inputfile);

% ----------------------------------------------------------------------------------------------------------
% Fundamental Constants 
q = 1.602e-19;                           % Charge on an electron [coul] 
eps0 = 8.85e-12;                         % Permittivity of vacuum [F/m] 
kB = 1.38e-23;                           % Boltzmann constant [J/K]
h = 6.63e-34;                            % Planck constant [J.s]
m0 = 9.11e-31;                           % Electron mass [kg]

% ----------------------------------------------------------------------------------------------------------
% Reading file etcetera
%%%inputfile = 'InputSLJ.txt';
fstring = fileread(inputfile);
fblocks = regexp(fstring,'\n','split');
Nb = length(fblocks);

% ----------------------------------------------------------------------------------------------------------
% Variable assignment
jj = 1;
vars = {};
vals = {};
for ii=1:Nb
   s = fblocks{ii};
   if(length(s)>1)
     s = s(find(~isspace(s)));   % Remove all spaces
     strings = split(s,'=');
     if(length(strings) > 1)
       vars{jj} = strings{1};
       vals{jj} = strings{2};
       jj = jj + 1;
     end
   end
end
Nv = length(vars);


% ----------------------------------------------------------------------------------------------------------
% Include additional inputfile (such as a template)
includefile = '';
for ii=1:Nv
  if(strfind(vars{ii},'&include'))
     includefile = vals{ii};
  end
end

if(length(includefile)>0);
  jj = 1;
  newfstring = fileread(includefile);
  newfblocks = regexp(newfstring,'\n','split');
  newNb = length(newfblocks);
  newvars = {};
  newvals = {};
  for ii=1:newNb
     s = newfblocks{ii};
     if(length(s)>1)
       s = s(find(~isspace(s)));   % Remove all spaces
       strings = split(s,'=');
       if(length(strings) > 1)
         newvars{jj} = strings{1};
         newvals{jj} = strings{2};
         jj = jj + 1;
       end
     end
  end
  vars = cat(2,newvars,vars);
  vals = cat(2,newvals,vals);
end
Nv = length(vars);



% ----------------------------------------------------------------------------------------------------------
% Default parameter values, then replaced below as needed
[sys,sc,inter] = DefaultInputVals();
names = split(inputfile,'.');
sys.RunDataFile = cat(2,names{1},'.out');


% ----------------------------------------------------------------------------------------------------------
% Parsing system (sys) parameters to find the temperature first
for ii=1:Nv
  % ......................................................................
  if(strfind(vars{ii},'sys.T'))
    sys.T = eval(vals{ii});
    Vt = sys.T*kB/q;
    sys.dV = Vt/2;   % Default bias step size
  end
end
Vt = sys.T*kB/q;

% ----------------------------------------------------------------------------------------------------------
% Parsing system (sys) parameters via if statements below
for ii=1:Nv

  % ......................................................................
  if(strfind(vars{ii},'sys.I0'))
    sys.I0 = eval(vals{ii});
  end
   
  % ......................................................................
  if(strfind(vars{ii},'sys.Vbi'))
    sys.Vbi = eval(vals{ii});
    Vt = sys.T*kB/q;
    if(sys.Vbi < Vt)
      sys.Vbi = Vt;
    end
  end

  % ......................................................................
  if(strfind(vars{ii},'sys.Vmin'))
    sys.Vmin = eval(vals{ii});
    Vt = sys.T*kB/q;
    if(sys.Vmin < (Vt/2))
      sys.Vmin = Vt/2;
    end
  end

  % ......................................................................
  if(strfind(vars{ii},'sys.Vmax'))
    sys.Vmax = eval(vals{ii});
    Vt = sys.T*kB/q;
    if(sys.Vmax < (Vt))
      sys.Vmax = Vt;
    end
  end


  % ......................................................................
  if((length(strfind(vars{ii},'sys.beta'))>0) && (length(vars{ii})<10))
    sys.beta = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sys.beta_phi'))
    sys.beta_phi = eval(vals{ii});
  end


  % ......................................................................
  if(strfind(vars{ii},'sys.Ngum'))
    sys.Ngum = eval(vals{ii});
    if(sys.Ngum < 50)  % Minimum bound 
      sys.Ngum = 50;
    end
  end

  % ......................................................................
  if(strfind(vars{ii},'sys.ct'))
    sys.ct = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sys.Nphi'))
    sys.Nphi = eval(vals{ii});
  end


end   % conclude system variable assignment loop

for ii=1:Nv   % reassign dV if specifically set in the input file
  if(strfind(vars{ii},'sys.dV'))
    sys.dV = eval(vals{ii});
    Vt = sys.T*kB/q;
    if(sys.dV > Vt) % Maximum voltage step for stability
      sys.dV = Vt; 
    end
  end
end

% ----------------------------------------------------------------------------------------------------------
% Parsing semiconductor (sc) parameters via if statements below
for ii=1:Nv

  % ......................................................................
  if(strfind(vars{ii},'sc.Nc'))
    sc.Nc = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.Nv'))
    sc.Nv = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.Eg'))
    sc.Eg = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.epsr'))
    sc.epsr = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.alpha'))
    sc.alpha = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.ND'))
    sc.ND = eval(vals{ii}); 
    if(sc.ND > (1e20*(100^3)));
      sc.ND = 1e20*(100^3);   % Maximum doping limit
    end
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.tau'))
    sc.tau = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.ratio_n'))
    sc.ratio_n = eval(vals{ii});
    if(sc.ratio_n>0.99)
      sc.ratio_n = 0.99;
    end
    if(sc.ratio_n<0.5)
      sc.ratio_n = 0.5;
    end     
  end

  % ......................................................................
  if(strfind(vars{ii},'sc.carrier_max'))
    sc.carrier_max = eval(vals{ii});
  end
end  % conclude semiconductor assignment loop

sc.ni = sqrt(sc.Nv*sc.Nc)*exp(-sc.Eg/(2*Vt));     % Intrinsic carrier concentration [#/m^3]
sc.eps = sc.epsr*eps0;                            % Overall dielectric constant
sc.p_o = (sc.ni^2)/sc.ND;                         % Bulk hole density [#/m^3]
sc.B = 1/(sc.ND*sc.tau);                          % Bimoleular coefficient benchamarked to n-type bulk following Gartner [m^3/s]
sc.Lp = sqrt( ((2*sc.eps)/q)*(1/sc.ND)*sc.Eg );   % Hole diffusion length [m]


% ----------------------------------------------------------------------------------------------------------
% Parsing interface (inter) parameters via if statements below
for ii=1:Nv

  % ......................................................................
  if(strfind(vars{ii},'inter.kp0'))
    inter.kp0 = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'inter.Lint'))
    inter.Lint = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'inter.Wp'))
    inter.Wp = eval(vals{ii});
  end

  % ......................................................................
  if(strfind(vars{ii},'inter.Lscreening'))
    inter.Lscreening = eval(vals{ii});
    if(inter.Lscreening < 5e-9)
      inter.Lscreening = 5e-9;
    end
  end

end  % conclude interface assignment loop

inter.Lsc = sqrt( ((2*sc.eps)/q)*(1/sc.ND)*sc.Eg ); % Semiconductor screening length when band bending is equal to the band gap 


% ========================================================================================================== %




