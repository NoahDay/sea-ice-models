if (exist('Qbutton')>0);                      
  scaleQ=1.0;                                 
  set(Qbutton,'String',scaleQ);               
end                                           
if (exist('Abutton')>0);                      
  A=203.3;                                    
  set(Abutton,'String',A);                    
end                                           
if (exist('Bbutton')>0);                      
  B=2.09;                                     
  set(Bbutton,'String',B);                    
end                                           
disp('okay');                                 
if (exist('Dbutton')>0);                      
  Dmag=0.44;                                  
  set(Dbutton,'String',Dmag);                
end                                           
disp('okay');                                 
if (exist('Tbutton')>0);                      
  Toffset=0.;                                 
  set(Tbutton,'String',Toffset);              
end                                           
if (exist('hadleybutton')>0);                 
  set(hadleybutton,'Value',0);                
  hadleyflag=0;                               
end                                           
if (exist('albedobutton')>0);                 
  set(albedobutton,'Value',0);                
  albedoflag=0;                               
end                                           
 if (exist('coldstartbutton')>0);             
  set(coldstartbutton,'Value',0);             
  coldstartflag=0;                            
end                                           
                                              
                                              