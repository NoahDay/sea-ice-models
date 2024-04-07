% Code up of North and Coakley's seasonal EBM model
% Simplified to eliminate the ocean domain
% Designed to run without a seasonal cycle and hence
% to stop once an equilibrium solution is reached.
% The model uses an implicit trapezoidal method
% so the timestep can be long.

%size of domain.
jmx=151;

% Choose parameters.
% scaleQ
if (exist('scaleQ')==0); scaleQ=1.; end

% OLR constant.
if (exist('A')==0); A=203.3; end

% OLR coef.
if (exist('B')==0); B=2.09; end

%heat diffusion coefficient.
if (exist('Dmag')==0); Dmag = 0.44; end

%heat diffusion coefficient.
Toffset=0.;
if (exist('coldstartflag')==1); 
  if (coldstartflag==1), Toffset = -40; end
end
  
%Simulate Hadley Cell with Lindzen and Farrell plan
if (exist('hadleyflag')==0); hadleyflag = 0.; end

%Remove albedo feedback
if (exist('albedoflag')==0); albedoflag = 0.; end

%heat capacity over land.
Cl = 0.2; % something small to make it equilibriate quickly

%time step in fraction of year
delt=1./50;
NMAX=1000; 

%set up x array.
delx = 2.0/jmx;
x = [-1.0+delx/2:delx:1.0-delx/2]';
phi = asin(x)*180/pi;

%obtain annual array of daily averaged-insolation.
%[insol] = sun(x);
%Legendre polynomial realizatin of mean annual insol.
Q = 338.5;
S = Q*(1-0.241*(3*x.^2-1)); 
S=scaleQ*S; S=S(:);

%set up inital T profile 
T = 20*(1-2*x.^2);
T=T(:);
T=T+Toffset;
%load T_final
Tinit=T;

%setup D(x) if simulating the Hadley Cell
%and calculate the matrix Mh and invM.
if (hadleyflag)
  xmp=[-1:delx:1];
  D=Dmag*(1+9*exp(-(xmp/sin(25*pi/180)).^6));
  D=D(:);
  [invM,Mh]=setupfastM(delx,jmx,D,B,Cl,delt);
else
  D=Dmag*ones(jmx+1,1);
  [invM,Mh]=setupfastM(delx,jmx,D,B,Cl,delt);
end

%Boundary conditions
%Set up initial value for h.
 alb=albedo(T,jmx,x,albedoflag);
 src   = (1-alb).*S/Cl-A/Cl; src=src(:);
 h=Mh*T+src;

%Global mean temperature
Tglob=mean(T);

% Timestepping loop
for n=1:NMAX
   Tglob_prev = Tglob;
    
% Calculate src for this loop.
   alb=albedo(T,jmx,x,albedoflag);
   src=((1-alb).*S-A)/Cl; src=src(:);

% Calculate new T.
   T=-invM*(0.5*(h+src)+T/delt);

% Calculate h for next loop.
   h=Mh*T+src;

% Check to see if global mean temperature has converged
   Tglob=mean(T);
   Tchange = Tglob-Tglob_prev;
   if (abs(Tchange) < 1.0e-12), break; end
end

%save T_final.mat T

% compute meridional heat flux and its convergence
a=6.37e+6; % earth radius in meters
[invM,Mh]=setupfastM(delx,jmx,D,0.,1.0,delt);
Dmp=0.5*( D(2:jmx+1)+D(1:jmx) );
divF=Mh*T;
F=-2*pi*a^2*sqrt(1-x.^2).*Dmp.*gradient(T,delx);


figure; 
subplot(3,1,1);
plot(phi,T,'.-','linewidth',1.5)
ylabel('Temperature'); xlabel('latitude');
set(gca,'position',[0.1300    0.71    0.7750    0.21]);
title(['Global mean temperature is ',num2str(Tglob,'%7.2f')]);
grid on;

subplot(3,1,2);
plot(phi,F*1e-15,'.-','linewidth',1.5)
ylabel('Poleward Heat Flux (10^{15} W)'); xlabel('latitude');
set(gca,'position',[0.1300    0.41    0.7750    0.21]);
grid on;

subplot(3,1,3);
plot(phi,divF,'.-',phi,(1-alb).*S,'o',phi,A+B*T,'.','linewidth',1.5)
ylabel('Energy Balance Terms (W m^{-2})'); xlabel('latitude');
set(gca,'position',[0.1300    0.130    0.7750    0.21]);
legend('\nabla F','SWd','LWu',4);
grid on;
u=axis; pos=u(3)-0.4*(u(4)-u(3));
text(-90,pos,['D = ',num2str(Dmag,'%7.2f'),...
       ',    Q/Qo = ',num2str(scaleQ,'%7.3f'),...
       ',    A = ',num2str(A,'%7.1f'),...
       ',    B = ',num2str(B,'%7.1f'),...
       ',    Toffset = ',num2str(Toffset,'%7.1f')] );

