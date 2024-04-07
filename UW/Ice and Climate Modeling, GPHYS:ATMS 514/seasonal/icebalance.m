   k=zeros(jmx,1); Fnet=k;
   k(ice)=conduct./h(ice);
   r(2*ice)=k(ice)*Tfrz-rprimew(ice);
   dev=zeros(2*jmx,1);
   dev(2*ice)=-Cw_delt+k(ice);
   Mt=M+diag(dev);
   I=inv(Mt)*r;
   T=I;
   T(2*ice)=min(Tfrz,I(2*ice));
   L=T(1:2:2*jmx);
   W=T(2:2:2*jmx);
   I=I(2:2:2*jmx);
   Fnet(ice)=Diff_Op(ice,:)*I-rprimew(ice)-B*W(ice)...
       -nu_fw(ice).*(I(ice)-L(ice));Fnet=Fnet(:);

% compute a heat flux from the ocean 

% the flux asymtotes to 2 W/m2 as the ice area 
% decreases down to a single gridbox of area (delx)
% and asymtotes to 0 as the ice area covers the hemisphere
% the ocean temperature is adjusted to 

nhice=ice(find(ice>jmx/2+1));         % nh grid cells with ice
shice=ice(find(ice<jmx/2));
nhocn=notice(find(notice>jmx/2+1));   % nh grid cells with ocn
shocn=notice(find(notice<jmx/2));
nhicearea=sum(fw(nhice));             % nh ice area
shicearea=sum(fw(shice));
nhmax=sum(fw(jmx/2+1:jmx));               % nh max possible ocn/ice area
shmax=sum(fw(1:jmx/2));
nhfw=2*min(2-2*(nhicearea-delx)/nhmax,2);   % nh fw under ice
shfw=2*min(2-2*(shicearea-delx)/shmax,2);
Fnet(nhice)=Fnet(nhice)+nhfw;
Fnet(shice)=Fnet(shice)+shfw;

nhocnarea=nhmax-nhicearea;                % nh ice-free ocn area
shocnarea=shmax-shicearea;
nhdW=nhfw*nhicearea./nhocnarea/Cw_delt;   % nh adjust water temp 
shdW=shfw*shicearea./shocnarea/Cw_delt;
W(nhocn)=W(nhocn)-nhdW;
W(shocn)=W(shocn)-shdW;
   