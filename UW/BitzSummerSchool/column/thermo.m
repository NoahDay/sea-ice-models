function [heat_added]=thermo(heat_added);
%
% heat budget
%
% need to read forcing Tair for entire grid
% add back in gairx,y
% deal with long dep of solar
%
% flux is positive towards the atm/ice or atm/ocean surface even if
% the flux is below the surface
%
% compute fluxes, temp and thickness over open water and ice
%

global centi
global rflsno rcpice  rcpsno  
global tsmelt tfrez 

global n1 nday dtau idter iyear iday
global hsstar
global fw hice hsnow ts tice eice esnow saltz
global hiout hsout tsout errout
global fsh flo upsens upltnt mualbedo io_surf snofal
persistent fid_20 fid_21 fid_23 fid_24 fid_30 fid_35

io=0;     % solar transmitted through top surface of ice
io1=0;    % solar absorbed by first layer
ib=0;     % solar transmitted through bottom surface of ice
fneti=0;   % net flux into top surface of ice from atmosphere and ice
fx=0;      % adjusted Fw accounting for times when ice melts away
fneg=0;    % net flux exchanged between ice and ocean
condb=0;   % conductive flux in ice at bottom surface
condt=0;   % conductive flux in ice at top surface
dq1=0;     % enthalpy change of first layer

delhib=0;        % ice thickness change at bottom for each category
delhit=0;        % ice thickness change at top for each category (excluding sublimation)
delhs=0;         % snow thickness change for each category (excluding sublimation)
subi=0;          % ice thickness change for each category from sublimation
subs=0;          % snow thickness change for each category from sublimation

layer=0;

% heat and energy accounting variables
e_init=0;
e_end=0;
heat_init=0;
heat_end=0;
difference=0;

%-----------------------------------------------------------------------;

% compute the initial enthalpy stored in ice./snow and mixed layer;
% and the heat added to check each step;
e_init=sumall;
heat_init=heat_added; 
fneg=0;

%fprintf(1,'%0.15g \n',e_init);

if (hice<5), 
  disp('model can not run without ice');
else
%-----------------------------------------------------------------------
% initialize ice temperature
%-----------------------------------------------------------------------
  tice(2:(n1+1))=gettmp;        % gettmp is a function

% for layer=1:3, fprintf(1,'%8.4f ',tice(layer+1)); end; fprintf(1,' A \n');
% eice=energ(tice(2:(n1+1)),saltz(2:(n1+1))); tice(2:(n1+1))=gettmp;
% for layer=1:3, fprintf(1,'%8.4f ',tice(layer+1)); end; fprintf(1,' B \n');


  if (hsnow<hsstar | esnow>0),  % wipe out small amount of snow
    fneg=snownrg;               % snownrg is a function
    hsnow = 0;
    tice(1)= tsmelt;
  end
  fneg=fneg/dtau; % turn fneg into a flux

  albedo=calc_albedo;
  fsh_net=fsh*(1-albedo); % net solar into ice surface from above
  fracsnow = hsnow/(hsnow + 0.1*centi); % patchy snow param
  io=0.;  % solar penetrating surf
  if (hsnow < hsstar), io=fsh_net*io_surf; end 

% fprintf(1,'%s %0.15g \n','A temp file',ts);
% for layer=0:n1, fprintf(1,'%0.15g ',tice(layer+1)); end; fprintf(1,'\n');
% fprintf(1,'%0.15g %0.15g %0.15g %0.15g %0.15g %0.15g \n',flo,io,fsh_net,upltnt,upsens,heat_added);

%-----------------------------------------------------------------------
% routine heart
%-----------------------------------------------------------------------
  [heat_added,fneti,condb,dq1,io1,ib,condt,ulwr]=tstmnew(flo,io,fsh_net,upltnt,upsens,heat_added);

% fprintf(1,'%s %0.15g \n','B temp file',ts);
% for layer=0:n1, fprintf(1,'%0.15g ',tice(layer+1)); end; fprintf(1,'\n');

 [delhib,delhs,delhit,subi,subs,fx]= growb(fneti,0.0,condb);

   hice=hice+delhit+delhib+subi;
   hsnow=hsnow+delhs+subs;
   fneg=fneg+fx;

%-----------------------------------------------------------------------
% update snow
%-----------------------------------------------------------------------
    if (snofal > 0.0)
      hs_init = hsnow;
      hsnow= max(hsstar,hsnow+snofal);
      dhs = hsnow-hs_init;
      tice(1) = (tice(1)*hs_init + tsmelt*dhs) / hsnow;
      heat_added = heat_added - dhs*rflsno;
    end

%-----------------------------------------------------------------------
% Energy budget diagnostics
%-----------------------------------------------------------------------
    heat_added=heat_added+fw*dtau;

    esnow= snownrg;      % after snofall added, snownrg is a function
    e_end= sumall;
    heat_end=heat_added;
    difference=((e_end-e_init)-(heat_end-heat_init))*0.001/dtau;



%-----------------------------------------------------------------------
% Output
%-----------------------------------------------------------------------
if (idter==nday)   % output some stuff

  nout=(iyear-1)*365+iday;
  hiout(nout)=hice;
  hsout(nout)=hsnow;
  tsout(nout)=ts;
  errout(nout)=difference;

%  if isempty(fid_20),
%   fid_20=fopen('outputs/thick','w');
%   fid_21=fopen('outputs/temp','w');
%   fid_23=fopen('outputs/fluxes','w');
%   fid_24=fopen('outputs/error','w');
%  end

  %error
%  fprintf(fid_24,'%0.15g %0.15g %0.15g \n',difference,...
%      (heat_end-heat_init)/1.e+9,(e_end-e_init)/1e+9);

  %fluxes
%  fprintf(fid_23,'%10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f \n',fsh_net,upltnt,upsens,ulwr,io,ib,io1,dq1,condt,condb);
 
  %temp
% fprintf(fid_21,'%8.5f ',ts);
% tice(2:(n1+1))=gettmp;
% for layer=0:n1, fprintf(fid_21,'%8.4f ',tice(layer+1)); end; fprintf(fid_21,'\n');

%% for layer=0:4, fprintf(fid_21,'%8.4f ',tice(layer+1)); end; fprintf(fid_21,'\n');
%% for layer=1:3, fprintf(fid_21,'%8.4f ',eice(layer)); end; fprintf(fid_21,'\n');

 %thick
% fprintf(fid_20,'%8.2f %8.2f %6.3f \n',hice,hsnow,albedo);

end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;


function [tmp]=gettmp;

% compute the midpoint temperature from the specific enthalpy for each layer;

global alpha gamma 
global rflice rcpice
global eice saltz
global n1

 layers=1:n1;
 q=eice(layers)+rflice-rcpice.*alpha.*saltz(layers+1);
 b=-q./rcpice;
 c=-gamma.*saltz(layers+1)./rcpice;

 b_2=b./2.0;
 tmp=-b_2-sqrt(b_2.*b_2-c);

end

