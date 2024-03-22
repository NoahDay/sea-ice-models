function [nrg]=energ(tmp,sal);
%
global rflice  rcpice
global alpha gamma 

% compute the specific enthalpy for non-new ice;
% relative to melting(negative quantity);
% assuming tmp is in celcius;

nrg = -rflice-rcpice.*(-alpha.*sal-tmp)-gamma.*sal./tmp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
