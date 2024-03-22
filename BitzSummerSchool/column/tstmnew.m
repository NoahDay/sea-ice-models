function [heat_added,f,condb,dq1,io1,ib,cond,ulwr]=tstmnew(flo,io,dswr,ultnt,usens,heat_added);
%
global n1 dtau hice hsnow  ts tice tw
global rcpice  rcpsno  
global sigma  esice
global tiny  hsmin  hsstar
global tffresh tmelt tsmelt tfrez qsmelt qmelt
global alpha gamma 
global kappa
global saltz


% This function calculates the evolution of the ice interior and surface
% temperature from the heat equation and surface energy balance;
%
% It solves the heat equation which is non-linear for saline ice;
% using an iterative newt-raps multi-equation solution.;
% scheme is either backwards or crank-nicholson, giving a tridiagonal;
% set of equations implicit in ti;
% Default differencing is backwards (CN seems to work, but the model;
% does not conserve heat when melting);
%
% must have a sufficient amount of snow to solve heat equation in snow;
% hsmin is the minimum depth of snow in order to solve for tice(0+1);
% if snow thickness < hsmin then do not change tice(0+1);
%
% the number of equations that must be solved by the tridiagonal solver;
% depends on whether the surface is melting and whether there is snow.;
% four soln_types are possible(see variable soln_type below):;
% 1 = freezing w./ snow, 2 = freezing w./ no snow,;
% 3 = melting w./ snow, and 4 = melting w./ no snow;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The call list I/O
% 
% vars sent
% tw;                  % water temp below ice
% flo;                 % long wave down at surface
% io; ib;              % solar penetrating top and bottom surfaces
% dswr;                % above surface net downward shortwave, indep of ts
% ultnt; usens;        % net upward latent, sensible
%                      % in this implementation with MU forcing, these are 
%                      % independent of temperature 

% vars sent and returned
% heat_added;          % diagnostic to keep track of work done on ice from atm.

% vars returned
f=0;                   % net flux at top surface including conductive flux in ice/snow
cond=0;                % conductive flux at top surface
condb=0;               % conductive flux at bottom surface
dq1=0; io1=0;          % solar absorbed in top layer and enthalpy change of top
ulwr = 0;              % net upward longwave flux

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local variables
ts_old=0; ti_old=zeros(1,n1+1);   % as above, except initially
dti=zeros(1,n1+2);   % incremental changes to ts and ti from newt-rhaps
ts_kelv=0;             % surface temperature in kelvin
z=zeros(1,n1+1);       % vertical coordinate of layer interfaces
qice=0;                % vapor at ice/snow surface 
ultnt_init=0; usens_init=0; ulwr_init=0;   % as above at start, depend on ts_old
ultnt_melt=0; usens_melt=0; ulwr_melt=0;   % as above at start, depend ts=melting
isol=zeros(1,n1+1);  % solar at layer interfaces
iabs=zeros(1,n1);    % solar absorbed in each layer
absorb=0;              % sum of Iabs
fo=0;                  % net flux at top surface excluding conductive flux in ice/snow
f_init=0; fo_init=0;   % as above, but init=depending on ts_old
fo_melt=0;             % as above, but melt=depending on ts=melting
dfo_dt=0;              % derivative of Fo wrt temperature
fofix=0;               % terms in Fo that are independent of temp
iru=0;                 % dummies used to compute F and Fo
condfix=0;             % terms in cond that are independent of temp
cond_init=0; condb_init=0;   % as above, but init=depending on ts_old
cond_melt=0;                 % as above, but melt=depending on ts=melting
specialk=0;            % modified conductivity for top layer of ice or snow
melts=0;               % the surface melting temp (either TSMELT or TMELT)
alph=0; bet=0;         % parameters for maintaining 2nd order accurate diff at boundar
% a,b,c are vectors that describe the diagonal and off-diagonal elements of;
% the matrix [a], such that [a] ti = r;
a=zeros(1,n1+2); b=zeros(1,n1+2); c=zeros(1,n1+2); d=zeros(1,1+2); 
r=zeros(1,n1+2);
ki=zeros(1,n1+1+1);  % layer conductivity divided by layer thickness
zeta=zeros(1,n1+1);  % the terms in heat equation that are independent of ti
delta=zeros(1,n1);   % gamma * Salinity / fresh ice heat capacity
eta=zeros(1,n1+1);   % time step / ice layer thickness / fresh ice heat capacity
cpi=zeros(1,n1);     % time step / ice layer thickness / fresh ice heat capacity
dh=0;                  % ice layer thickness
dt_dh=0; dt_hs=0;      % time step / ice or snow thickness
n=0;                   % number of equations solved by tridiag solver
n1p2=n1+2;
layers=1:n1;           % counters for ice layers
layers0=0:n1;          % counters for snow and ice layers
iloop=0;                 % counter for iterations of numerical loop
iter=0;                % counter for number of times try to find a numerical soln
keepiterating=true;    % flag indicating newton-rhaps. gone astray
errit=0;               % the absolute value of the maximum dti
errmax=1.0e-6 ;        % the maximum error allowed for dti in degrees
errmax=1.0e-4 ;        % the maximum error allowed for dti in degrees

theta=1.0;             % 1.0 for backwards and 0.5 for CN (DO NOT LET THETA=0.)
soln_type=0;      % indicates the case of snow vs. no snow and melting vs. freezing


% -------------------------------------------------------------------------;
% setup helpful parameters;
% -------------------------------------------------------------------------;
dh=hice./n1;
dt_dh=dtau./dh;
dt_hs=0.0; if(hsnow>hsmin), dt_hs=dtau./hsnow; end;
ts_old=ts;  ti_old=tice;
tbot=min(tw,tmelt);
[ki]=conductiv(tice,tbot,hsnow,dh,n1);
%for layer=0:n1, fprintf(1,'%0.15g ',ki(layer+1)); end; fprintf(1,'\n');

% the solar radiation absorbed internally;
z=cumsum([0 dh*ones(1,n1)]);
isol=exp(-kappa.*z);
iabs=io.*(isol(layers)-isol(1+layers));

if(hsnow>hsmin),
  alph = 2.0.*(2..*hsnow + dh) ./(hsnow+dh);
  bet = -2.0.*hsnow.*hsnow./(2..*hsnow+dh)./(hsnow+dh);
  condfix=ki(0+1).*(alph.*tice(0+1)+bet.*tice(1+1));
  specialk=ki(0+1).*(alph+bet);
  melts=tsmelt;
  ulwr_melt = esice.*(tsmelt+tffresh)^4-flo;
  qice=qsmelt;
else;
  alph = 3.;
  bet = -1../3.;
  specialk=ki(1+1).*8.0./3.0;
  condfix=ki(1+1).*(alph.*tice(1+1)+bet.*tice(2+1));
  melts=tmelt;
  ulwr_melt = esice.*(tsmelt+tffresh)^4-flo;
  qice=qmelt;
end;

% -------------------------------------------------------------------------;
% get a first guess for ts based on ti from previous step;
% assume the surface temperature is at the melting point and;
% see if it is balanced. note that this alters ts_old;
% as far as the cn method goes;
% this could be debated - though it doesn't affect energy conservation.;
% -------------------------------------------------------------------------;
cond_melt=condfix-specialk.*melts;
fofix=dswr-io+flo-ultnt-usens;
fo_melt=dswr-io-ulwr_melt-ultnt-usens;
f =fo_melt+cond_melt;   % assuming ts=melting
if(f<0.0) ;
  f=0.;
  ts=iter_ts_mu(specialk,fofix,condfix);
else;
  ulwr_init =ulwr_melt;
end;

condb_init=ki(n1+1+1).*( 3..*(tbot-tice(n1+1)) -(tbot-tice(n1-1+1))./3.0 );
f_init=f;
ts_old=ts;

% -------------------------------------------------------------------------;
% beginning of iterative proceedure;
% -------------------------------------------------------------------------;
iter=1;  iloop=1; errit = 10.; keepiterating = true;

while (keepiterating),

while ((errit>errmax) & (iloop<20)),   % numerical loop
  % -------------------------------------------------------------------------;
  % setup terms that depend on the cpi and;
  % initial temp(ti_old and ts_old);
  % -------------------------------------------------------------------------;
  cpi(layers)=rcpice+gamma.*saltz(layers+1)./tice(layers+1)./ti_old(layers+1);
  eta(0+1)=0.0;
  if(hsnow>hsmin), eta(0+1)=dtau./(hsnow.*rcpsno);  end;
  eta(layers+1)=dt_dh./cpi(layers);
  
  if(theta<1) 
  else;
    zeta=ti_old;
    zeta(1+layers)=  zeta(1+layers) + eta(1+layers).*iabs;
  end;

  f =fo_melt+cond_melt;    % F assuming ts=melting
  if(f<=0),     % solve heat equation for ice./snow using flux bc;
    ts_kelv=ts+tffresh;
    iru= esice.*(ts_kelv).^4;
    dfo_dt=-4.*esice.*ts_kelv.^3;
    fo = fofix-iru;
    if(hsnow>hsmin) ;
      % -------------------------------------------------------------------------;
       soln_type=1; % case of freezing with snow layer
      % -------------------------------------------------------------------------;
      [a,b,c,r]=getabc(tice,tbot,zeta,delta,ki,eta,n1,1);
      cond=ki(0+1).*(alph.*(tice(0+1)-ts)+bet.*(tice(1+1)-ts));
      a(0+2)=-eta(0+1).*ki(0+1).*(alph+bet);
      c(0+2)=eta(0+1).*(bet.*ki(0+1)-ki(1+1));
      b(0+2)=1+eta(0+1).*(alph.*ki(0+1)+ki(1+1));
      r(0+2)=-zeta(0+1)+a(0+2).*ts+b(0+2).*tice(0+1)+c(0+2).*tice(1+1);
      a(-1+2)=0.;
      c(-1+2)=-ki(0+1).*alph;
      d(-1+2)=-ki(0+1).*bet;
      b(-1+2)=-dfo_dt-c(-1+2)-d(-1+2);
      r(-1+2)=-fo-cond;
      % this is the row operation to get rid of d(-1+2);
      b(-1+2)=c(0+2).*b(-1+2)-d(-1+2).*a(0+2);
      c(-1+2)=c(0+2).*c(-1+2)-d(-1+2).*b(0+2);
      r(-1+2)=c(0+2).*r(-1+2)-d(-1+2).*r(0+2);
      n=fix(n1+2);
      [dti]=tridag(a,b,c,r,n);
      ts=ts-dti(-1+2);
      tice(layers0+1)=tice(layers0+1)-dti(layers0+2);
  
    else;
    % -------------------------------------------------------------------------;
      soln_type=2;  % case of freezing with no snow layer
    % -------------------------------------------------------------------------;
      [a,b,c,r]=getabc(tice,tbot,zeta,delta,ki,eta,n1,2);
      a(1+2)=-eta(1+1).*ki(1+1).*(alph+bet);
      c(1+2)=-eta(1+1).*(ki(2+1)-bet.*ki(1+1));
      b(1+2)=1+eta(1+1).*(ki(2+1)+alph.*ki(1+1));
      r(1+2)=-zeta(1+1)+a(1+2).*ts+b(1+2).*tice(1+1)+c(1+2).*tice(2+1);
      a(0+2)=0.;
      c(0+2)=-ki(1+1).*alph;
      d(0+2)=-ki(1+1).*bet;
      b(0+2)=-dfo_dt-c(0+2)-d(0+2);
      r(0+2)=-fo-ki(1+1).*(alph.*(tice(1+1)-ts)+bet.*(tice(2+1)-ts));
      % this is the row operation to get rid of d(0+2);
      b(0+2)=c(1+2).*b(0+2)-d(0+2).*a(1+2);
      c(0+2)=c(1+2).*c(0+2)-d(0+2).*b(1+2);
      r(0+2)=c(1+2).*r(0+2)-d(0+2).*r(1+2);
      n=fix(n1+1);
      [dti(2:n1p2)]=tridag(a(2:n1p2),b(2:n1p2),c(2:n1p2),r(2:n1p2),n);
      ts=ts-dti(0+2);
      tice(layers+1)=tice(layers+1)-dti(layers+2);

    end
  else
  % these dont need to iterated but it is okay to do so
    ts=melts;
    if(hsnow>hsmin) ;
  % -------------------------------------------------------------------------;
      soln_type=3;   % case of melting with snow layer
  % -------------------------------------------------------------------------;
      [a,b,c,r]=getabc(tice,tbot,zeta,delta,ki,eta,n1,1);
      a(0+2)=-eta(0+1).*ki(0+1).*(alph+bet);
      c(0+2)=eta(0+1).*(bet.*ki(0+1)-ki(1+1));
      b(0+2)=1+eta(0+1).*(alph.*ki(0+1)+ki(1+1));
      r(0+2)=-zeta(0+1)+a(0+2).*ts+b(0+2).*tice(0+1)+c(0+2).*tice(1+1);
      a(0+2)=0.0;

%      for layer=2:n1p2, fprintf(1,'%0.15g ',a(layer)); end; fprintf(1,'\n');
%      for layer=2:n1p2, fprintf(1,'%0.15g ',b(layer)); end; fprintf(1,'\n');
%      for layer=2:n1p2, fprintf(1,'%0.15g ',c(layer)); end; fprintf(1,'\n');
%      for layer=2:n1p2, fprintf(1,'%0.15g ',r(layer)); end; fprintf(1,'\n');

      n=fix(n1+1);
      [dti(2:n1p2)]=tridag(a(2:n1p2),b(2:n1p2),c(2:n1p2),r(2:n1p2),n);
      tice(layers0+1)=tice(layers0+1)-dti(layers0+2);
%      fprintf(1,'iter iloop %0.15g %0.15g \n', iter,iloop); fprintf(1,'%0.15g ',ts);
%      for layer=0:n1, fprintf(1,'%0.15g ',tice(layer+1)); end; fprintf(1,'\n');
%      keyboard;
    else;
  % -------------------------------------------------------------------------;
      soln_type=4;   % case of melting with no snow layer
  % -------------------------------------------------------------------------;
      [a,b,c,r]=getabc(tice,tbot,zeta,delta,ki,eta,n1,2);
      a(1+2)=-eta(1+1).*ki(1+1).*(alph+bet);
      c(1+2)=-eta(1+1).*(ki(2+1)-bet.*ki(1+1));
      b(1+2)=1+eta(1+1).*(ki(2+1)+alph.*ki(1+1));
      r(1+2)=-zeta(1+1)+a(1+2).*ts+b(1+2).*tice(1+1)+c(1+2).*tice(2+1);
      a(1+2)=0.0;
      n=fix(n1);
      [dti(3:n1p2)]=tridag(a(3:n1p2),b(3:n1p2),c(3:n1p2),r(3:n1p2),n);
      tice(layers+1)=tice(layers+1)-dti(layers+2);
    end;
  end;
  
  % -------------------------------------------------------------------------;
  % end numerical proceedure, see if need to reiterate;
  % -------------------------------------------------------------------------;
  
  errit=0.0;
  if(hsnow>hsmin),errit=abs(dti(0+2)); end;
  errit=max([errit,abs(dti(layers+2))]);
  ts=min(ts,melts);
  if(hsnow>hsmin) ;
    condfix=ki(0+1).*(alph.*tice(0+1)+bet.*tice(1+1));
    cond_melt=condfix-specialk.*tsmelt;
  else;
    condfix=ki(1+1).*(alph.*tice(1+1)+bet.*tice(2+1));
    cond_melt=condfix-specialk.*tmelt;
  end;
  
%  fprintf(1,'tstm errit %0.15g \n', errit); fprintf(1,'%0.15g ',ts);
%  for layer=0:n1, fprintf(1,'%0.15g ',tice(layer+1)); end; fprintf(1,'\n');
  
  iloop=fix(iloop+1);

end % numerical loop

% -------------------------------------------------------------------------;
% check to see if ti < melting, if not try;
% numerical loop again with different initial conditions;
% -------------------------------------------------------------------------;
keepiterating=false;

if any(  tice(0+1)>(tsmelt+tiny) | tice(layers+1) >-alpha.*saltz(layers+1) ), 
  keepiterating=true; 

  if (iter==1) ;
    % this works 99 times out of 100
    ts=ts_old-20.;
    tice(layers+1)=ti_old(layers+1)-20.;
  elseif (iter==2) ;
    % this works the rest of the time
    ts=melts-0.5;
    tice(layers+1)=tmelt-0.5;
  else;
    % when this error occurs, the model does not conserve energy;
    % it should still run, but there probably is something wrong;
    % that is causing this sceme to fail;
    keepiterating=false; % give up on finding a solution
    fprintf(1,'%s \n','WARNING converges to ti>TMELT');
    fprintf(1,'%0.15g ',ts_old);
    for l=(1):(n1), fprintf(1,'%0.15g ',ti_old(l+1)); end;
    fprintf(1,'%0.15g \n',tbot);
    fprintf(1,'%0.15g ',h);
    fprintf(1,'%0.15g ',dswr);
    fprintf(1,'%0.15g ',flo);
    fprintf(1,'%0.15g \n',io);
    ts=melts;
    tice(layer0+1)=min(-alpha.*saltz(layer0+1),tice(layer0+1));
  end;

  if(hsnow>hsmin) ;
    condfix=ki(0+1).*(alph.*tice(0+1)+bet.*tice(1+1));
    cond_melt=condfix-specialk.*tsmelt;
  else;
    condfix=ki(1+1).*(alph.*tice(1+1)+bet.*tice(2+1));
    cond_melt=condfix-specialk.*tmelt;
  end;
end;

iter=fix(iter+1);
iloop=0; % reset numerical loop counter

end;


% -------------------------------------------------------------------------;
% continue from here when done iterative
% finish up by updating the surface fluxes, etc.;
% -------------------------------------------------------------------------;
if(ts<melts) ;
  cond=(condfix-specialk.*ts);
  fo=-cond;
  f=0.0;
  ulwr = esice.*(ts+tffresh).^4-flo;
  if(theta<1) ;
    ulwr = theta.*ulwr +(1.-theta).*ulwr_init;
  else;
    % gaurantees the surface is in balance;
    ulwr = cond+dswr-io-ultnt-usens;
  end;
else;
  fo=fo_melt;
  f =fo_melt+cond_melt;
  cond=cond_melt;
  if(theta<1) ;
    ulwr = theta.*ulwr_melt +(1.-theta).*ulwr_init;
  else;
    ulwr = ulwr_melt;
  end;
end;

if(errit>errmax) ;
% when this error occurs, the model does not conserve energy;
% it should still run, but there probably is something wrong;
% that is causing this sceme to fail;
fprintf(1,'%s ','WARNING No CONVERGENCE');
fprintf(1,'%0.15g ',errit);
fprintf(1,'%0.15g ',i);
fprintf(1,'%0.15g \n',j);
fprintf(1,'%0.15g ',ts,tice);
fprintf(1,'%0.15g ',ts_old,ti_old);
fprintf(1,'%0.15g \n',tbot);
fprintf(1,'%0.15g ',h);
fprintf(1,'%0.15g ',dswr);
fprintf(1,'%0.15g ',flo);
fprintf(1,'%0.15g \n',io);
  ts=tmelt; tice(layers0+1)=tmelt;
  f=0.0; fo=0.0; ultnt = 0.0;
end;

% condb is positive if heat flows towards ice;
condb=ki(n1+1+1).*( 3..*(tbot-tice(n1+1)) -(tbot-tice(n1-1+1))./3.0 );
absorb=sum(iabs(1:n1));
ib=io-absorb;
% solar passing through the bottom of the ice
condb=theta.*condb+(1.-theta).*condb_init;
f=theta.*f+(1.-theta).*f_init;
fo=theta.*fo+(1.-theta).*fo_init;
heat_added=heat_added+(fo+absorb)*dtau;
dq1=(gamma.*saltz(1+1).*(1./tice(1+1)-1./ti_old(1+1))+...
     rcpice.*(tice(1+1)-ti_old(1+1)))./dt_dh;
io1=iabs(1);
if(hsnow>hsmin) ;
  cond=(ki(1+1+1).*(tice(1+1+1)-tice(1+1))-ki(1+1).*(tice(1+1)-tice(1-1+1)));
else;
  cond=(ki(1+1+1).*(tice(1+1+1)-tice(1+1)))-cond;
end;

end

% -------------------------------------------------------------------------;

function [a,b,c,r]=getabc(ti,tbot,zeta,delta,k,eta,ni,lfirst);

alph=3.0; bet=-1.0./3.0 ;
% if there is snow lfirst=1 otherwise it is 2;
layers=lfirst:(ni-1);
a(layers+2)=-eta(layers+1).*k(layers+1);
c(layers+2)=-eta(layers+1).*k(layers+1+1);
b(layers+2)=1-c(layers+2)-a(layers+2);
r(layers+2)=-zeta(layers+1)+a(layers+2).*ti(layers-1+1) ...
  +b(layers+2).*ti(layers+1)+c(layers+2).*ti(layers+1+1);
a(ni+2)=-eta(ni+1).*(k(ni+1)-bet.*k(ni+1+1));
c(ni+2)=0.0;
b(ni+2)=1+eta(ni+1).*(k(ni+1)+alph.*k(ni+1+1));
r(ni+2)=-zeta(ni+1)+a(ni+2).*ti(ni-1+1)+b(ni+2).*ti(ni+1) ...
 -eta(ni+1).*(alph+bet).*k(ni+1+1).*tbot;
end

% -------------------------------------------------------------------------;

function [u]=tridag(a,b,c,r,n);
global n1

bet=0; gam=zeros(1,n1+2);
% if(b(1)==0.) pause 'tridag: rewrite equations';
bet=b(1);  u(1)=r(1)./bet;

for layer=2:n;
  gam(layer)=c(layer-1)./bet;
  bet=b(layer)-a(layer).*gam(layer);
  % if(bet==0.0) pause 'tridag failed';
  u(layer)=(r(layer)-a(layer).*u(layer-1))./bet;
end;

for layer=n-1:-1:1;
  u(layer)=u(layer)-gam(layer+1).*u(layer+1);
end; 
end

% -------------------------------------------------------------------------;

function [ki]=conductiv(ti,tbot,hsnow,dh,ni);
global kimin beta kappai kappas 
global saltz tice hsmin

tmax=-0.1 ;

layer=2:ni;
ki(layer+1)=kappai+beta.*(saltz(layer+1)+saltz(layer+1+1)) ...
    ./min(tmax,(tice(layer-1+1)+tice(layer+1)));

ki(ni+1+1)=kappai+beta.*saltz(ni+1+1)./tbot;
ki(1+1)=kappai+beta.*saltz(1+1)./min(tmax,tice(1+1));

ki=max(ki,kimin);
ki=ki./dh;

ki(0+1)=0.0;
if(hsnow>hsmin) ;
  ki(0+1)=kappas./hsnow;
  ki(1+1)=2.0.*ki(1+1).*ki(0+1)./(ki(1+1)+ki(0+1));
end;

end

% -------------------------------------------------------------------------;

function [ts]=iter_ts_mu(k,fofix,condfix);
global tffresh esice tiny

% first guess
ts=-20;
iter=0;
keepiterating=true;

while (keepiterating),
  ts_kelv=ts+tffresh;
  iru= esice.*ts_kelv.^4;
  cond=condfix-k.*ts;
  df_dt=-k-4.*esice.*ts_kelv.^3;
  f = fofix-iru+cond;
  if(abs(df_dt)<tiny),
    keepiterating=false;
  else
    dt=-f./df_dt;
    ts=ts+dt;
    iter=fix(iter+1);
    if(iter>20 | (abs(dt)>0.001)), keepiterating=false; end;
  end
end

end




