qqq=VEH_v*3.6/1.609;
   for i=1:101 
        if VEH_v(i)>=60*1.609/3.6
        t_60=i;
        break;
        end
   end
   t_60=t_60-3+(60-qqq(t_60-1))/(qqq(t_60)-qqq(t_60-1))
   
   
    plot(t, VEH_v*3.6/1.609, 'r');    % Actual velocity
    xlabel('Time [s]')
    ylabel('Velocity [MPH]')
    title('velocity vs time ')
    grid
    
    
%     
%    for i=1:101
%       if VEH_x(i)>=0.25*1609
%           t_qm=i;
%            v_qm=VEH_v(i);
%           break; 
%       end
%    end
% t_qm=t_qm-2+(0.25*1609-VEH_x(t_qm-1))/(VEH_x(t_qm)-VEH_x(t_qm-1))
% v_qm=VEH_v(i-1)+(0.25*1609-VEH_x(i-1))/(VEH_x(i)-VEH_x(i-1))*(VEH_v(i)-VEH_v(i-1));
% v_qm=v_qm*3.6/1.609
% 
% % 4-4 0-20 20-40 30-50 50-70 mph time intervals
% t1=0;
% t2=0;
% t3=0;
% t4=0;
% t5=0;
% for i=1:101
%    if t1==0
%        if VEH_v(i)>=20*1.609/3.6
%            t1=i;
%        end 
%    end
%    if (t1>0&&t2==0)
%        if VEH_v(i)>=30*1.609/3.6
%            t2=i;
%        end
%    end
%    if (t1>0&&t2>0&&t3==0)
%        if VEH_v(i)>=40*1.609/3.6
%            t3=i;
%        end
%    end
%    if (t1>0&&t2>0&&t3>0&&t4==0)
%        if VEH_v(i)>=50*1.609/3.6
%            t4=i;
%        end
%    end
%    if (t1>0&&t2>0&&t3>0&&t4>0&&t5==0)
%        if VEH_v(i)>=70*1.609/3.6
%            t5=i;
%            break;
%        end
%    end
% end
% 
% t1=t1-2+(20-qqq(t1-1))/(qqq(t1)-qqq(t1-1));
% t2=t2-2+(30-qqq(t2-1))/(qqq(t2)-qqq(t2-1));
% t3=t3-2+(40-qqq(t3-1))/(qqq(t3)-qqq(t3-1));
% t4=t4-2+(50-qqq(t4-1))/(qqq(t4)-qqq(t4-1));
% t5=t5-2+(70-qqq(t5-1))/(qqq(t5)-qqq(t5-1));
% t1
% t3-t1
% t4-t2
% t5-t4