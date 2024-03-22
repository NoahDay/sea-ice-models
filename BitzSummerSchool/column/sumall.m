function [energall]=sumall;
%
global  hice hsnow eice esnow n1

energall=sum(eice)*hice/n1 + esnow * hsnow;

end

