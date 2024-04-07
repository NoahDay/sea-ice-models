global guihnd

if exist('guihnd'), set(guihnd,'handlevisibility','off'); end

plotno=3;
pltfig=figure; 

axhnd=axes('parent',pltfig,'position',[0.13 0.7012 0.3270 0.2238],...
 'XLim',[0 1],'Ylim',[-90 90]);
plhnd=plot(fl,phi,'parent',axhnd); 
title(['Geography for ',casename,' Case'],'parent',axhnd)
text(0.05,-20,'LAND','rotation',90);text(0.8,-20,'OCEAN','rotation',90,'parent',axhnd);
xlabel('land/ocean fraction','parent',axhnd)
ylabel('latitude','parent',axhnd);

axhnd=axes('parent',pltfig,'position',[0.578 0.7012 0.3270 0.2238],...
	   'XLim',[0 360],'Ylim',[-90 90]); gca=axhnd;
contourf(1:360,phi,insol); shading flat;
title('Insolation','parent',axhnd);
hold on; [c,h] = contour(1:360,phi,insol,'w-');clabel(c,h); hold off; 
xlabel('Day of Year','parent',axhnd); ylabel('latitude','parent',axhnd);

days=(nstepinyear-1):-1:0; tt=n-days; 
cmin=min(min(min(L_out(:,tt))),min(min(W_out(:,tt))));
cmax=max(max(max(L_out(:,tt))),max(max(W_out(:,tt))));
cmin=floor(cmin/5)*5;
clevs=cmin:5:cmax;
axhnd=axes('parent',pltfig,'position',[0.1300    0.4056    0.3270    0.2238]);
gca=axhnd; contourf(3:6:360,phi,L_out(:,tt),clevs); 
shading flat;title('Land Temperature','parent',axhnd);
hold on; [c,h] = contour(3:6:360,phi,L_out(:,tt),clevs,'w-');clabel(c,h); hold off; 
xlabel('Day of Year','parent',axhnd); ylabel('latitude','parent',axhnd);

axhnd=axes('parent',pltfig,'position',[0.5780    0.4056    0.3270    0.2238]);
gca=axhnd; contourf(3:6:360,phi,W_out(:,tt),clevs); 
shading flat;title('Ocean Temperature','parent',axhnd);
hold on; [c,h] = contour(3:6:360,phi,W_out(:,tt),clevs,'w-');clabel(c,h); hold off; 
xlabel('Day of Year','parent',axhnd); ylabel('latitude','parent',axhnd);

axhnd=axes('parent',pltfig,'position',[0.1300    0.1100    0.3270    0.2238]);
gca=axhnd;
if (ice_model), 
  y=mean(h_out(:,tt)');
  miny=min(h_out(:,tt)');
  maxy=max(h_out(:,tt)');
  fill([phi; flipud(phi)],[miny fliplr(maxy)],'c')  
  ylabel('Sea Ice Thickness','parent',axhnd);
  xlabel('latitude','parent',axhnd);
  u=axis; axis([-90 90 0 u(4)]);
else
 axis([-90 90 0 1]);
 text(-50,0.8,'No explicit sea ice model','parent',axhnd) 
end

u=axis; pos=u(3)-0.25*(u(4)-u(3));
text(-90,pos,['D = ',num2str(mean(Dmag),'%7.2f'),...
       ',    Q/Qo = ',num2str(scaleQ,'%7.3f'),...
       ',    A = ',num2str(A,'%7.1f'),...
       ',    B = ',num2str(B,'%7.1f'),...
       ',    ecc = ',num2str(ecc,'%7.4f'),...
       ',    obl = ',num2str(obl,'%7.2f'),...
       ',    per = ',num2str(per,'%7.1f')] ,'parent',axhnd);
  
if (ice_model), 
  icy=zeros(size(h_out));j=find(h_out>0.1);  icy(j)=1;
  % area of each latitude circle is 1/jmx of globe
  earthrad=6.37e6;
  gridarea=4*pi*earthrad^2/jmx;
  j=find(phi>0);  extentNH=gridarea*sum(icy(j,:));
  j=find(phi<=0); extentSH=gridarea*sum(icy(j,:));
  axhnd=axes('parent',pltfig,'position',[0.578 0.11 0.327 0.1072]);
  plot(tday/360,extentNH,tday/360,extentSH,'--','parent',axhnd)
  ylabel('Sea Ice Extent','parent',axhnd);xlabel('years','parent',axhnd)
  u=axis;
  text(0.2*u(2),u(3)+0.75*(u(4)-u(3)),'NH solid, SH dashed');
  axhnd=axes('parent',pltfig,'position',[0.578 0.2516 0.327 0.1072]);
else
  axhnd=axes('parent',pltfig,'position',[0.578 0.11 0.327 0.2238]);
end
Tg=(L_out'*fl+W_out'*fw)/jmx;
days=((1/delt)-1):-1:0; tt=n-days;
plot(tday/360,Tg,'parent',axhnd);ylabel('Global Mean Temperature','parent',axhnd);
if not(ice_model), xlabel('years','parent',axhnd); end
u=axis;
text(0.6*u(2),u(3)+0.05*(u(4)-u(3)),['T_g = ',num2str(mean(Tg(tt)))],'parent',axhnd);

orient tall;

if exist('guihnd'), set(guihnd,'handlevisibility','on'); end
