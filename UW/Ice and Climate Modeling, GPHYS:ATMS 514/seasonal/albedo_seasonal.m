function [alb_l,alb_w]=albedo2(L,W,x);

% recalculate albedo.


%k=find(W<-2);  alb_w(k)=0.68;
%k=find(L<-10); alb_l(k)=0.7;



alb_w=0.313+0.08*(3*x.^2-1)/2-.05;
alb_l=0.313+0.08*(3*x.^2-1)/2+.05;

k=find(W<=-2);  alb_w(k)=0.6;
k=find(L<=-2); alb_l(k)=0.6;

