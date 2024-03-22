function init_column
%
global centi
global sigma  esice
global rflice  rflsno  rslice  rslsno  rcpice  rcpsno  
global tffresh tmelt tsmelt tfrez qsmelt qmelt
global alpha gamma 
global kappa kimin beta kappai kappas 

global n1 nyrs
global tiny  hsmin  hsstar
global frzpt fw area hice hsnow tw tbot ts tice eice esnow saltz
global hiout hsout tsout errout

sigma=5.67e-5;      % Stefan constant
esice=sigma;        % ice emissivity times Stefan constant

rhoice = 917.e-03;  % density of ice
rhosno = 330.e-03;  % density of snow

cpice = 2054.e+04;  % heat capacity of fresh ice
cpsno=2113.e+04;    % heat capacity of snow

vlocn= 2.501e10;    % latent heat of vaporization freshwater
slice= 2.835e10;    % latent heat of sublimation freshwater
flice=slice-vlocn;  % latent heat of fusion freshwater

rflice=flice*rhoice; % specific latent heat of fushion ice
rflsno=flice*rhosno; % specific latent heat of fushion snow
rslice=slice*rhoice; % specific latent heat of sublim ice
rslsno=slice*rhosno; % specific latent heat of sublim snow
rcpice=cpice*rhoice; % specific heat capacity of fresh ice
rcpsno=cpsno*rhosno; % specific heat capacity of snow


tffresh = 273.16;  % Freezing Temperature Of Freshwater
tmelt   = -0.0;
tsmelt  = -0.00;   % Melting Temperature Of The Snow
tfrez = 271.2-tffresh; % Approx Freezing Temperature Of The Ocean

% coefficients for computing saturation vapor pressure
ai=21.8746;       % over ice
bi=7.66-273.16;   % over ice
qs1=0.622*6.11/1013.; %(mol weight of water:dry air)/(surface pressure in mb)*6.11
qsmelt=qs1*exp(ai*tsmelt/(tsmelt-bi));
qmelt =qs1*exp(ai*tmelt/(tmelt-bi));

% parameter for computing melting temp as a function of salinity.
% i.e., melting=Tffresh-alpha*salinity
alpha=0.054;

kappa = 1.5e-2; % Solar Extinction Coef In Ice
kimin = 0.1e+5; % Minimum Conductivity In Ice
beta = 0.1172e+5; % Conductivity, K=Kappai+Beta*Salinity/T
gamma = rflice*alpha; % Heat Capacity C=Cpi+Gamma*Salinity/T**2

% for thin ice growth in the lead, want to use predetermined
% values for kappai, based on estimate for -5 deg C surface temp
kappai=2.034e+5; % thermal conductivity of fresh ice
kappas=0.31e+5; % thermal conductivity of snow

tiny=1e-6;

% if there is not at least a centimeter of snow do not bother
% with conductive effects of snow in the heat equation
% do not rub out though, give it a chance to accumulate
% also allow it to melt
% albedo depends on snow depth so small amounts of snow hanging around
% will not screw-up the albedo
% WHEN THE SNOW IS <HSMIN ITS TEMPERATURE IS EQUAL TO THE SURFACE
hsmin=1.;  hsstar=0.0005;


% first guess at initial hice and area and set misc. variables to 0;
% ice is uniformly set to 2.8 m with no open water;
centi=100.0;kilo=1000.;
frzpt=tfrez;
fw=2.0.*kilo;
area=1.0;
hice=2.53*centi;
hsnow=0.2827*centi;
tw=frzpt;
ts=tfrez;

% variables for timeseries 
hiout=zeros(1,nyrs*365);
hsout=zeros(1,nyrs*365);
tsout=zeros(1,nyrs*365);
errout=zeros(1,nyrs*365);

saltz=salinity_prof(n1);

if (n1==10) ;
  ts =-29.2894;
  tice(0+1)=-23.1128;
  tice(1+1)=-15.7925;
  tice(2+1)=-14.1042;
  tice(3+1)=-12.4534;
  tice(4+1)=-10.8423;
  tice(5+1)=-9.2777;
  tice(6+1)=-7.7687;
  tice(7+1)=-6.3255;
  tice(8+1)=-4.9605;
  tice(9+1)=-3.6879;
  tice(10+1)=-2.5230;
else;
  tice(0+1)=-23.16-(tfrez+23.16)./n1;
  for layer=1:n1
    tice(layer+1)=-23.16+(layer-1).*(tfrez+23.16)./n1;
  end;
end;
tbot=min(tw,tmelt);

esnow= snownrg;
layers=1:n1;
eice = energ(tice(layers+1),saltz(layers+1));

end

