% code up of North and Coakley's seasonal EBM model
% with sea ice
%clear;
seasonal_setup;

r=zeros(2*jmx,1);
if (ice_model),
%set up initial sea ice thickness of 2 m where freezing
  ice=find( W<Tfrz);
  notice=find( W>=Tfrz);
  h=zeros(jmx,1); k=h;
  h(ice)=2;
  if (albedoflag),
    alb_l=clim_alb_l(:,thedays(1)); alb_w=clim_alb_w(:,thedays(1));
  else
    [alb_l,alb_w]=albedo_seasonal(L,W,xfull);
  end    
  S=insol(:,ts);
  rprimew=A-(1-alb_w).*S;
  rprimel=A-(1-alb_l).*S;
  r(1:2:2*jmx,1)=L*Cl_delt-rprimel;
  r(2:2:2*jmx,1)=W*Cw_delt-rprimew;
  icebalance;
else
  ice=[];
end

L_out = zeros(jmx,length(n_out));
W_out = L_out; 
if( ice_model), h_out=W_out; end

clear tday yr day
% Begin loop
idx_out=1;
for n = 1:nts
%for n = 1:10
   tday(n) = ts +2 +1+(n-1)*360*delt;
   yr(n)=floor((-1+tday(n))/360);
   day(n)=floor( tday(n)- yr(n)*360);
   
%create initial albedo.
  if (albedoflag),
    nn = day(n);
    nn=nn-360*delt; if (nn<0), nn=360-360*delt *0.5; end
    alb_l=clim_alb_l(:,nn); alb_w=clim_alb_w(:,nn);
  else
    [alb_l,alb_w]=albedo_seasonal(L,W,xfull);
  end    

% calculate insolation
   S = insol(:,day(n));

%Source terms
  rprimew=A-(1-alb_w).*S;
  rprimel=A-(1-alb_l).*S;
  r(1:2:2*jmx,1)=L*Cl_delt-rprimel;
  r(2:2:2*jmx,1)=W*Cw_delt-rprimew;

 if ( ice_model  ),
% first consider where sea ice already exists
  ice=find(h>0.001);
  notice=find(h<=0.001);
  if (length(ice)>0), 
   icebalance;
   h(ice)=h(ice)-delt_Lf*Fnet(ice);
   h(ice)=max(0.,h(ice));
  end

% second consider the conditions of new ice growth over the ocean  
  cold=find( T(2*notice)<Tfrz);
  new=notice(cold);
  if (length(new)>0),
    h(new)=-Cw/Lfice*( W(new)-Tfrz );
    W(new)=Tfrz;
  end

 else
   T=invM*r;
   L=T(1:2:2*jmx);
   W=T(2:2:2*jmx);
 end

%output 
   if n==n_out(idx_out) 
     L_out(:,idx_out) = L(:);
     W_out(:,idx_out) = W(:);
     if( ice_model), h_out(:,idx_out) = h(:); end
     idx_out=idx_out+1;
   end
end
n=idx_out-1;
days=(nstepinyear-1):-1:0; tt=n-days; 
thedays=day(tt);
Lann=L_out(:,tt);
Wann=W_out(:,tt);
%save temperatures.mat  Lann Wann thedays

plotoutput
