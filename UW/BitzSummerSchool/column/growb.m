function [delb,delhs,delh,subi,subs,alarm]=dh(fneti,ultnt,condb);
%
global rflsno rcpsno
global tsmelt

global n1 nday dtau
global tiny  hsmin  hsstar
global frzpt fw hice hsnow tbot ts tice eice esnow saltz


alarm=false;
alarm2=false;

esnow=snownrg;
eice=energ(tice(2:(n1+1)),saltz(2:(n1+1)));
enet=sumall;

dhi=hice/n1*ones(1,n1);
dhs=hsnow;
delh=0.; delhs=0.; delb=0.;
subi=0; subs=0;

if (ultnt > 0),
 % be warned that this portion of the matlab port
 % has not been tested
  etop=ultnt*dtau;
  evnet=enet+etop-rvlsno*hsice-rvlice*hice;
  if (evnet>0),
    subi=-h; subs=-hs;
    alarm = true;
  else
    [subi,subs,alarm2]=surfsub(enet,etop,dhs,dhi,delh,delhs);
    dhs=hs+subs;
    si=subi;
    for layer=1:n1
      s=max(-dhi(layer),si);
      dhi(layer)=dhi(layer)+s;
      si=si-s;
    end
    if alarm2, alarm=true; disp('subsidence error'); end
  end
end

if (fneti>0 & ~alarm),
% it is no big deal to melt at surf when Ts is below melting
% because it never melts more than a hundreth of an Angstom
  etop=fneti*dtau;
  enet=enet+etop;
  if (enet > 0.0)
    delh=-(hice+subi);
    delhs=-(hsice+subs);
    fx=condb-enet/dtau;
    alarm=true;
    disp(['melted though all layers from top']);
  else
    [delh,delhs,alarm2]= surfmelt(etop,dhs,dhi,delh,delhs);
    dhs=dhs+delhs;
    si=delh;
    for layer=1:n1
       s=max(-dhi(layer),si);
       dhi(layer)=dhi(layer)+s;
       si=si-s;
    end
    if alarm2, alarm=true; disp('surfmelt error'); end
  end
end

if (~alarm),
  delht=delh+subi;
  fx=fw;
  ebot=dtau*(fw-condb);
  if (ebot < 0)
    egrow=energ(tbot,saltz(n1+2));
    delb=ebot/egrow;
  else
    % melt at bottom, melting temperature = Tbot
    egrow=0.; % on purpose
    if ((enet+ebot)>0.)
      delb=-(hice+delht);
      delhs=-(hsnow+delhs+subs);
      fx=condb-enet/dtau;
      keyboard
      disp('melted though all layers from bottom');
      alarm=true;
    else
      [delb,delh,delhs,alarm2]= botmelt(ebot,dhs,dhi,delh,delhs,delb);
      if alarm2, disp('problem with botmelt'); alarm=true; end
    end
  end
  if ~alarm, eice=adjust(egrow,delb,delht); else, keyboard; end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [delh,delhs,alarm] = surfsub(enet,etop,dhs,dh,delh,delhs);
global rvlsno eice esnow n1

u=esnow-rvlsno;
finished=false; alarm=false;
if (u*dhs< 0.0) 
  % convert etop into equivalent snowmelt
  delhs = etop/u
  if ((dhs + delhs)<0.) 
     % Melt only some of the snow
     etop=0.0;
     enet=enet+esnow*delhs;
     finished=true;
  else
    % Melt all of the snow and some ice too
    delhs=-dhs;
    etop=etop+u*dhs;
    enet=enet+es*delhs;
  end
end

if ~finished
  for layer=1:n1
    u=eice(layer)-rvlice;
    if (-u*dh(layer) <= etop) ;
      % melt partial layer
      delh=delh+etop/u;
      enet=enet+etop/u*ei(layer);
      etop=0.0
      finished=true
    else, 
      % melt out whole layer
      delh=delh-dh(layer)
      etop=etop+u*dh(layer)
      enet=enet+dh(layer)*eice(layer)
    end
  end
  if ~finished, disp(['ERROR in surfsub',etop]); alarm=true; end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delh,delhs,alarm] = surfmelt(etop,dhs,dh,delh,delhs);
global eice esnow n1 tice tsmelt rflsno rcpsno

finished=false; alarm=false;
u=esnow;  
if (u*dhs<0.0)
  % convert etop into equivalent snowmelt
  delhs = etop/u;
%  keyboard
  if ((dhs + delhs)>0.) 
    % Melt only some of the snow
    etop=0.0;
    finished=true; 
  else
    % Melt all of the snow and some ice too
    delhs=-dhs;
    etop=etop+u*dhs;
  end
end

if ~finished
  for layer=1:n1
    u=eice(layer);
    if (-u*dh(layer)>=etop) 
       delh=delh+etop/u;
       etop=0.0;
       finished=true; 
    else
       delh=delh-dh(layer);
       etop=etop+u*dh(layer);
    end
  end
  if ~finished, disp(['ERROR in surfmelt',etop]); alarm=true; end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delb,delh,delhs,alarm] = botmelt(ebot,dhs,dh,delh,delhs,delb);
global eice esnow n1 hsnow

finished=false; alarm=false;

%for layer=1:n1, fprintf(1,'%8.4f ',eice(layer)); end; fprintf(1,'\n');
%fprintf(1,'ebot %0.15g %0.15g \n',ebot,delb);
%for layer=1:n1, fprintf(1,'%8.4f ',dh(layer)); end; fprintf(1,'\n');

for layer=n1:-1:1
  u=eice(layer);
  if (-u*dh(layer)>=ebot) 
    delb=delb+ebot/u;
    ebot=0.0;
    finished=true; 
  else
    delb=delb-dh(layer);
    ebot=ebot+u*dh(layer);
  end
end
%keyboard

if ~finished,
  % finally melt snow if necessary
  u=esnow;
  if (-u*dhs>=ebot) 
    delhs=delhs+ebot/u;
    ebot=0.0;
    finished=true;
  else
    disp('melted completely through all ice and snow');
    disp(['ERROR in botmelt']); alarm=true;
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function [e_tw] = adjust(egrow,delb,delh);

% * Adjusts temperature profile after melting/growing.
% *
% * eice is the energy density after updating tice from
% * the heat equation without regard delh and delb
% *
% * hice is the thickness from previous time step!
% * h_tw is the NEW thickness 
% *
% * delb is negative if there is melt at the bottom
% * delh is negative if there is melt at the top
% *
% * generally _tw is a suffix to label the new layer spacing variables
% *
global tiny hice n1 eice

e_tw=eice;
if ~((abs(delb)<tiny) & (delh>-tiny))

  h_tw=hice+delb+delh;

  if (h_tw <= 0.0), 
    e_tw=zeros(1,n1); 
  else

    % layer thickness
    delta=hice/n1;  delta_tw=h_tw/n1;

    % z is positive down and zero is relative to the top
    % of the ice from the old time step
    z=zeros(1,n1+2); z_tw=zeros(1,n1+1);
    z_tw(1)=-delh;
    layers=2:n1;
    z(layers)=delta*(layers-1);
    z_tw(layers)=z_tw(1)+delta_tw*(layers-1);

    z(n1+1)=hice;
    z(n1+2)=hice+max(delb,0.0);
    z_tw(n1+1)=z_tw(1)+h_tw;

    fract=zeros(n1,n1+1);
    for l_tw=1:n1
      for  l=1:(n1+1)
        fract(l_tw,l)=( min(z_tw(l_tw+1),z(l+1))- ...
             max(z_tw(l_tw),z(l)) );
      end
    end
    fract=fract/delta_tw;   fract=max(fract,0.0);

    e_tw=[eice egrow]*fract';

  end
end
end
