function nrg=snownrg
%
% the energy required to heat the snow to;
% the melting temperature and melt it;
% currently dont heat./cool the meltwater to the lead temperature;
% the snow layer melts at tsmelt;

global rflsno  rcpsno
global hsnow tice
global tsmelt

nrg=-rflsno+rcpsno.*(tice(1)-tsmelt);

end
