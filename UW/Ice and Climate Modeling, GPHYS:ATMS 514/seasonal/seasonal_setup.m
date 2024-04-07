% DEFAULTS
defaults

if (exist('jmx')==0); jmx=jmxdef; end
if (exist('runlength')==0); runlength=runlengthdef; end
if (exist('scaleQ')==0); scaleQ=scaleQdef; end
if (exist('A')==0); A=Adef; end
if (exist('B')==0); B=Bdef; end
if (exist('Dmag')==0); Dmag = Dmagdef; end
if (exist('nu')==0); nu = nudef; end
if (exist('Cl')==0); Cl = Cldef; end
if (exist('Cw')==0); Cw = Cwdef; end
if (exist('coldstartflag')==0), coldstartflag = coldstartdef; end
if (exist('hadleyflag')==0); hadleyflag = hadleyflagdef; end
if (exist('albedoflag')==0); albedoflag = albedoflagdef; end
if (exist('obl')==0);  obl  = obldef;      end
if (exist('ecc')==0);  ecc  = eccdef; end
if (exist('per')==0); per = perdef;      end
if (exist('ice_model')==0);  ice_model = ice_modeldef; end
if (exist('land')==0);  land = landdef; end
if (exist('casename')==0); casename=casenamedef; end

%size of domain.
jmx=2*floor(jmx/2);

% freezing temperature of ocean in deg C
% 83 is ratio of Lf/Cp and 50 is depth of ML
Tfrz=-2;
conduct=2;
Lfice=9.8*83.5/50; % this is about right

%time loop parameters.
ts=90; % first day of run is present day equinox
tf=runlength-0.25;
nstepinyear=60;
delt=1./nstepinyear;
nts=floor(tf/delt);
n_out=1:nts;

%set up x array.
delx = 2.0/jmx;
x = [-1.0+delx:delx:1.0-delx]';
xfull=(-1+delx/2:delx:1-delx/2)';
phi = asin(xfull)*180/pi;

%set up inital T profile
if coldstartflag, Toffset=-40; else, Toffset = 0.; end
L = 7.5+20*(1-2*xfull.^2)+Toffset;
W = 7.5+20*(1-2*xfull.^2)+Toffset;

%heat diffusion coefficient.
if (hadleyflag)
  D=Dmag*(1+9*exp(-(x/sin(25*pi/180)).^6));
else
  D=Dmag*ones(length(x),1);
end

%land fraction.
fl = 0.05*ones(jmx,1);
if sum(land(1:5)=='Preca')==5,
  % pangea (middle Cambrian)
  j=find(phi<=-45); fl(j)=0.3;
  j=find(phi>-45 & phi<=30); fl(j)=0.5;
elseif sum(land(1:5)=='Ordov')==5
  % gondwanaland
  j=find(phi<=-60); fl(j)=0.95;
  j=find(phi>-60 & phi<=70); fl(j)=0.3;
elseif sum(land(1:5)=='Symme')==5
  % Symmetric
  fl(1:jmx)=0.34;
else % modern
  fl = 0.38*ones(jmx,1);
  j=find(phi<=-60); fl(j)=0.95;
  j=find(phi>-60 & phi<=-40); fl(j)=0.05;
  j=find(phi>-40 & phi<=20); fl(j)=0.25;
  j=find(phi>20 & phi<=70); fl(j)=0.5;
end
fl=fl*0.34/mean(fl);

%ocean fraction.  
fw = 1-fl;


%obtain annual array of daily averaged-insolation.
[insol] = seasonal_solar(xfull,obl,ecc,per);
insol=scaleQ*insol(:,[360 1:359]);


Cw_delt=Cw/delt;
Cl_delt=Cl/delt;
delt_Lf=delt/Lfice;

nu_fw = nu./fw;
nu_fl = nu./fl;

lambda=D/delx/delx.*(1-x.*x);
a=[0; -lambda];
c=[-lambda; 0];
b=-a-c;
Diff_Op =-( diag(b) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1) ) ;

bw=Cw_delt+B+nu_fw-a-c;
bl=Cl_delt+B+nu_fl-a-c;

Mw = diag(bw) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1) ;
Ml = diag(bl) + diag(c(1:jmx-1),1) + diag(a(2:jmx),-1) ;

M=zeros(2*jmx,2*jmx);
for j=1:jmx
  M(2*j-1,1:2:2*jmx)=Ml(j,:);
  M(2*j-1,2*j)=-nu_fl(j);
  M(2*j,2:2:2*jmx)  =Mw(j,:);
  M(2*j,2*j-1)=-nu_fw(j);
end

invM=inv(M);

% make climatological albedo if no albedo feedback
if (albedoflag),
  load temperatures.mat 
  clim_alb_l=zeros(jmx,360);
  clim_alb_w=clim_alb_l;
  n=1;
  for t = thedays
    [clim_alb_l(:,t),clim_alb_w(:,t)]=albedo_seasonal(Lann(:,n),Wann(:,n),xfull);
    n=n+1;
  end
end
