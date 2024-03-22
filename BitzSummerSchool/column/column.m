function [hiout, hsout, tsout, errout]=driver(LWpert,Nyr,timeofyear);
% [HIOUT, HSOUT, TSOUT, ERROUT]=COLUMN(LWPERT,NYR,TIMEOFYEAR)  
%   Driver for Bitz and Lipscomb column model.
%   If no options are given the model runs with default 
%   values for inputs LWPERT, NYR, and TIMEOFYEAR.
%
%   Be sure to type a semicolon after calling this routine.
%   
%   For example, type: 
%   [hiout, hsout] = column(-2,50);    
%   or
%   [hiout, hsout] = column(2,20,2);
%   or
%   column(8);
%
%   Intended values for the options are 
%      LWPERT downward longwave perturbation:  -10 to 10 (default 
%          is zero), units are W/m2
%      NYRS   run length: 10 to 150 (default is 20), units are years
%      TIMEOFYEAR portion of year when LWPERT is applied, 0, 1, or 
%          2 = allyear, winter only, or summer only, respectively
%          (default is 0).  Winter only is from autumnal to vernal 
%          equinox.
%   The user may specify LWPERT alone, LWPERT and NYR alone, or all 
%   three options.
%  
%   The output variables are timeseries of daily values for
%      HIOUT ice thickness in cm
%      HSOUT snow thickness in cm
%      TSOUT surface temperature in deg C
%      ERROUT error in the energy conservation converted into W/m2
%
%   A figure is plotted with the daily timeseries
%  
%
%   C.M. Bitz, June 24, 2007

global n1 nyrs nday dtau idter iyear iday
global fsh flo upsens upltnt mualbedo io_surf snofal
global hiout hsout tsout errout
global centi

if (nargin==0), LWpert=0.; end;
if (nargin>1),  nyrs=Nyr; else, nyrs=20;  end;    % run length; 
if (nargin<3), timeofyear=0; end;  % flag for applying LWpert 0=allyear, 1=winter, 2=summer

n1=10;   % number of layers
nday = 2;
dtau = 86400/nday;    % time step, note that nday is increased in last 3 yrs of the run

init_column;  % initialize a bunch of other global vars

% declare vars for this routine
start_time=cputime;
firststep=true;  
fsh_n1=0; flo_n1=0; dnsens_n1=0; dnltnt_n1=0; mualbedo_n1=0;
fsh_n=0;  flo_n=0;  dnsens_n=0;  dnltnt_n=0;  mualbedo_n=0;
e_init=0; e_end=0;  % energy in the ice and snow
heat_added=0;       % running total of heat added to the ice and snow


% compute the initial energy in the ice./snow and mixed layer;
e_init=sumall;

% data from fletcher 1965, used by maykut and untersteiner;
% interpolated from monthly to daily without preserving monthly means. 
% This is less than ideal
fid=fopen('data.mu71','r');
data=reshape(fscanf(fid,'%g'),5,365)';
fclose(fid);

pertdays=ones(1,365); equin=172;
if (timeofyear==1);  pertdays(:,112:223)=0;  
elseif (timeofyear==2);   pertdays(:,[1:111 224:365])=0; end

for iyear=1:nyrs;
 if (iyear>nyrs-3),
  nday = 6;
  dtau = 86400/nday;    % time step
 end
  for iday=1:365;
%    fprintf(1,'%s %0.15g\n', 'day=',iday);

    % prepare to interpolate the forcing data
    % n1 = today, n = yesterday;
    fsh_n=fsh_n1;
    flo_n=flo_n1;
    dnsens_n=dnsens_n1;
    dnltnt_n=dnltnt_n1;
    mualbedo_n=mualbedo_n1;
    fsh_n1=data(iday,1); 
    flo_n1=data(iday,2); 
    dnsens_n1=data(iday,3);  
    dnltnt_n1=data(iday,4);
    mualbedo_n1=data(iday,5);

    if firststep,
      fsh_n=fsh_n1;
      flo_n=flo_n1;
      dnsens_n=dnsens_n1;
      dnltnt_n=dnltnt_n1;
      mualbedo_n=mualbedo_n1;
      firststep=false;
    end;

    for idter=1:nday;

    % linear spline daily forcing data forml timestep=day./nday;
      fsh=fsh_n+(fsh_n1-fsh_n).*idter./nday;
      flo=LWpert*1000*pertdays(iday)+flo_n+(flo_n1-flo_n).*idter./nday;
      upsens=dnsens_n+(dnsens_n1-dnsens_n).*idter./nday;
      upltnt=dnltnt_n+(dnltnt_n1-dnltnt_n).*idter./nday;
      mualbedo=mualbedo_n+(mualbedo_n1-mualbedo_n).*idter./nday;
      mualbedo=mualbedo-0.0475;

      io_surf=0.3;
      snofal=centi*snowfall(iday)/nday;
      heat_added=thermo(heat_added);

    end;
  end;
  fprintf(1,'%s %0.15g\n', 'finished year ',iyear);
end;


e_end=sumall;
end_time=cputime;

fprintf(1,'\n%s \n','Energy Totals from the Run converted into  W/m^2');
fprintf(1,'%s ','energy change of system =');
fprintf(1,'%0.15g \n',(e_end-e_init)*0.001/(nyrs*86400*365));
fprintf(1,'%s ','heat added to the ice/snow');
fprintf(1,'%0.15g \n',heat_added*0.001/(nyrs*86400*365));
fprintf(1,'%s ','-->  difference: ');
fprintf(1,'%0.15g \n',(e_end-e_init-heat_added)*0.001/(nyrs*86400*365));
fprintf(1,'%s %4.1f s\n\n','run time',end_time-start_time);

fprintf(1,'%s\n','Final Year Statistics');
htme=(nyrs-1)*365+(1:365); tme=(nyrs-1)*365+(32:91);
fprintf(1,'%s %6.2f\n','Mean Thickness',mean(hiout(htme)));
fprintf(1,'%s %7.3f\n','Mean Feb-Mar Temperature',mean(tsout(tme)));

clf; tme=(1:length(hiout))/365;
subplot(2,2,1);plot(tme,hiout); xlabel('year'); ylabel('ice thickness - cm');
subplot(2,2,2);plot(tme,hsout); xlabel('year'); ylabel('snow depth - cm');
subplot(2,2,3); plot(tme,tsout); xlabel('year'); ylabel('surface temperature - C');
subplot(2,2,4); plot(tme,errout); xlabel('year'); ylabel('error - W m^{-2}');

end


function snow=snowfall(id);
% snowfall from maykut and untersteiner 1971;
snow=0;
if(id<=119 | id>=302) ;
  snow= 2.79e-4;
elseif (id>=120 & id<=150) ;
  snow= 1.61e-3;
elseif(id>230) ;
  snow= 4.16e-3;
else;
  snow= 0.0;
end;

end
