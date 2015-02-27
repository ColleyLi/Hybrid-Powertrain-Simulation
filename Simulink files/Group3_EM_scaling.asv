
%================================================================================================
% EM scaling
   load('em_scale.mat')
   load('EM2_map')
   % pme, cm dimensionless model
   
   
   % input desired torque
   % EM maxT
    T_scale=100;
   %%%% scale
   rl_0=0.0338
   rl=rl_0*(T_scale)^(2/3)/(180)^(2/3)
   
   %rl=0.0338 % r x l
    
   % keep the geometry propostional
   k=0.0338/(0.3^2); %r=k*l, rl=k*l^2
   l=sqrt(rl/k);
   r=k*l;

   T_s=2*pi*(r^2)*l*map_np
   w_s=cm/r
   
  figure(1)
  C_eta = [ 0.7 0.8 0.85 0.875 0.89 0.90];
  [c1,h1] = contour(w_s',T_s',map_ns',C_eta );
  clabel(c1,h1);
  xlabel('Speed (rad/s)')
  ylabel('Torque (N.m)')
  title('Efficiency')
  hold on
  grid on 
  title('Motor/Alternator efficiency map ')
  axis([0 1000 -40 40])
  
  %save('em_scale.mat','map_np','cm','map_ns','maxT_w')
  
  T_0=48
  l_0=0.3
  
  max_T2(1:15)=T_0*(l^3)/(l_0^3)
  
    r_0=0.5*0.225
%    l=0.3
  
  maxT_w2=maxT_w*r_0/r
 
  p_peak=max_T2(15)*maxT_w2(15)
  

  
  for n=15:29
      max_T2(n)=p_peak/maxT_w2(n)
  end
  
  plot(maxT_w2,max_T2,'r')
  plot(maxT_w2,-max_T2,'r')
  
   maxTw_max=max(maxT_w2)
   xp = [0 maxT_w2' maxTw_max];
   yp = [700 max_T2 700];
   patch(xp,yp,'w')
   
   yp = [-700 -max_T2 -700];
   patch(xp,yp,'w')
   
   EM2.maxT_w=maxT_w2;
   EM2.maxT_T=max_T2';
   EM2.minT_T=-max_T2';
   EM2.map_w=w_s;
   EM2.map_T=T_s;
   EM2.map_eff=map_ns;  
   
   maxp=zeros(1,length(maxT_w));
   for i=1:length(maxT_w)
      maxp(i)=maxT_w(i)*max_T2(i);
   end
% % Maxp_EM=max(maxp);
% % EM2.m=(Maxp_EM/1000/105)^(3/2)*45;
% % EM2.map=Maxp_EM;
% %    save('EM2_map','EM2')
%================================================================================================
