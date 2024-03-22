function saltz=salinity_prof(n);
%
% ice salinity
saltmax=3.2; saltmin=1.0;
for layer=1:n;
  zrel=(layer-0.5)./n;
  saltz(layer+1)=saltmax./2..*...
    ( 1+sin(pi.*(zrel.^(0.40706205./(zrel+0.57265966))-0.5)) );
end;
saltz(0+1)=0.0;       % snow layer salinity, not used
saltz(n+2)=saltmax;   % base of the ice, probably not used

% Could try this to mimic Winton model
% if n==2, 
%  saltz(2)=3.2;
%  saltz(3)=0.;
% end

% fprintf(1,'%s \n','salt profile');
% for layer=1:n, fprintf(1,'%0.15g \n',saltz(layer+1)); end;
