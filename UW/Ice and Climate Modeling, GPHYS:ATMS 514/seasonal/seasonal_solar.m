function [insol] = sun(xi,obl,ecc,long)

% code for insolation from Wen-Wei Pan, based on Berger 
% formulas (JAS, 35, 1978)

%Some parameters
npts=length(xi);
%units of W/m2 and Daily averaged insolation
t1=2808/2.0754;
dr_conv=pi/180;
rd_conv=1/dr_conv;

%time array and extra stuff.
nts=360;
ti=[1:(360-1)/(nts-1):360];
distance = (1-ecc^2)./(1+ecc*cos(dr_conv*(450+ti-long)));
s_delt=-sin(dr_conv*obl)*cos(dr_conv*ti);
c_delt=sqrt(1.0-s_delt.^2);
t_delt = s_delt./c_delt;
delt=asin(s_delt)*rd_conv;

phi = asin(xi)*rd_conv;

wk = zeros(npts,nts);

%loop over latitudes and days to produce insolation
for i = 1:nts
   for j = 1:npts
      if delt(i) > 0.0                     
         if phi(j) >= 90-delt(i)
            wk(j,i)=t1*xi(j)*s_delt(i)/(distance(i)^2);
         elseif ((-phi(j) >= (90 - delt(i))) & (phi(j)<0))
            wk(j,i) = 0;
         else
	    c_h0 = -tan(dr_conv*phi(j))*t_delt(i);
            h0   = acos(c_h0);
            wk(j,i) = t1*(h0*xi(j)*s_delt(i)+cos(dr_conv*phi(j))...
                         *c_delt(i)*sin(h0))...
                          /(distance(i)^2*pi);
         end
      else
	 if phi(j) >= (90 +delt(i))
            wk(j,i)=0;
         elseif (-phi(j)>=(90+delt(i)) & (phi(j)<0))
            wk(j,i)=t1*xi(j)*s_delt(i)/distance(i)^2;
         else
            c_h0=-tan(dr_conv*phi(j))*t_delt(i);
            h0 = acos(c_h0);
            wk(j,i) = t1*(h0*xi(j)*s_delt(i)+cos(dr_conv*phi(j))...
	                 *c_delt(i)*sin(h0))...
                         /(pi*distance(i)^2);
         end
      end
   end
end


% now interpolate onto full resolution array
%insol=interp2(ti,xi,wk,[1:360],x);
insol=wk;
