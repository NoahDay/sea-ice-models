defaults

if (exist('jmxbutton')>0);                  
  set(jmxbutton,'String',jmxdef);                
  jmx=jmxdef;
end                                       
if (exist('runlengthbutton')>0);                  
  set(runlengthbutton,'String',runlengthdef);                
  runlength=runlengthdef;
end                                       
if (exist('Qbutton')>0);                  
  set(Qbutton,'String',scaleQdef); 
  scaleQ=scaleQdef;
end                                       
if (exist('Abutton')>0);                  
  set(Abutton,'String',Adef);                
  A=Adef;
end                                       
if (exist('Bbutton')>0);                  
  set(Bbutton,'String',Bdef);                
  B=Bdef;
end                                       
if (exist('Dbutton')>0);                  
  set(Dbutton,'String',Dmagdef);                
  Dmag=Dmagdef;
end                                       
if (exist('eccbutton')>0);                  
  set(eccbutton,'String',eccdef);                
  ecc=eccdef;
end                                       
if (exist('oblbutton')>0);                  
  set(oblbutton,'String',obldef);                
  obl=obldef;
end                                       
if (exist('perbutton')>0);                  
  set(perbutton,'String',perdef);                
  per=perdef;
end                                       
if (exist('nubutton')>0);                  
  set(nubutton,'String',nudef);                
  nu=nudef;
end                                       
if (exist('Clbutton')>0);                  
  set(Clbutton,'String',Cldef);                
  Cl=Cldef;
end                                       
if (exist('Cwbutton')>0);                  
  set(Cwbutton,'String',Cwdef);                
  Cw=Cwdef;
end                                       

if (exist('hadleybutton')>0);             
  set(hadleybutton,'Value',1);  
  hadleyflag=hadleyflagdef;
end                                       
if (exist('albedobutton')>0);             
  set(albedobutton,'Value',0);            
  albedoflag=albedoflagdef;                           
end                                       
if (exist('coldstartbutton')>0);         
  set(coldstartbutton,'Value',0);         
  coldstartflag=coldstartdef; 
end                                       
if (exist('icemodelbutton')>0);         
  set(icemodelbutton,'Value',0);         
  ice_model=ice_modeldef;                        
end                                       
if (exist('casenamebutton')>0);         
  casename=get(casenamebutton,'String');
  if length(casename)<1,
   set(casenamebutton,'String',casenamedef);
   casename=casenamedef;
  end
end                                       
                                          
  