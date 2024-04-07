function [albedo]=calc_albedo
%
global hsmin hsnow ts tmelt

pert= -0.0;
albedo=0.63+pert;

if(hsnow>hsmin); 
  albedo = 0.8+pert;
  if ( ts>= (tmelt-0.01) ), albedo=0.75+pert; end;
end;

end
