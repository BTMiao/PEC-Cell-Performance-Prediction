% ========================================================================================================== %
% DISCLAIMER                                                                                                 %
% The authors and publishers make no warranties about the software, and disclaim liability for all uses of   %
% the software, to the fullest extent permitted by applicable law.  The authors do not recommend use of this %
% software for any purpose.  It is made available, solely to clarify points made in their analysis of        %
% photoelectrochemical (PEC) devices.  When using or citing the software, you should not imply endorsement   %
% by the authors.                                                                                            %
%                                                                                                            %
% DESCRIPTION                                                                                                %
% o Performs the LU decomposition and provides the results of the  matrix inversion                          %
%      [A][x]=[fn] => [LU][x]=[fn]                                                                           %
%   Inputs: an (lower diagonal), bn (diagonal), cn (upper diagonal), fn                                      %
%                                                                                                            %
% ========================================================================================================== %

function [n] = LUdec(an,bn,cn,fn)

    % ------------------------------------------------------------------------------------------------------
    Np=length(bn);
    alphan(1) = bn(1);
    for i=1:Np-1
        betan(i)=an(i)/alphan(i);
        alphan(i+1)=bn(i+1)-betan(i)*cn(i);
    end

    % ------------------------------------------------------------------------------------------------------
    % Solution of Lv = f 
    vn(1) = fn(1);
    for i = 2:Np
        vn(i) = fn(i) - betan(i-1)*vn(i-1);
    end

    % ------------------------------------------------------------------------------------------------------
    % Solution of Ux = v 
    tempn = vn(Np)/alphan(Np);
    n(Np)=tempn;
    for i = (Np-1):-1:1
        tempn = (vn(i)-cn(i)*n(i+1))/alphan(i);
        n(i) = tempn;
    end

end  % Function End

% ========================================================================================================== %


