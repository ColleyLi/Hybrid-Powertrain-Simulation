
%============================================================================================
% EMSC Control Strategy
load EM2_map.mat
load Mat_files/TRAN_Data.mat
load Mat_files/EM_Data.mat;

Qhlv=43.0E06;    % gasoline Fuel Lower Heating Value [J/kg]
soc_h=0.85
soc_l=0.55

% size of look up table
m1=12 %x Treq
m2=10 %y  ice_w_e
m3=5  %z soc_e

em_trq_e=zeros(m2+1,m1+1,m3) % initalize matirx
% soc f penalty
%x_soc=(soc-(soc_h+soc_l)/2)/((soc_h-soc_l)/2)
x_soc=0;
x_soc0=linspace(-1,1,5);

%em data
new_em_maxT_w=EM2.maxT_w
new_em_maxT_T=EM2.maxT_T
new_em_minT_T=EM2.minT_T
new_em_map_T=EM2.map_T
new_em_map_w=EM2.map_w
new_em_map_eff=EM2.map_eff
rtl=-2; %maximum engine2EM charging torque

for gear=1:4 % select gear

 
 for n3=1:m3 % soc sweep
    x_soc=x_soc0(n3)
    weight_c=1;
    f_penalty=(1-(x_soc.^3)/2)*weight_c; % penalty
     
    
    for n2=0:m2 % speed sweep
        ice_w_e=n2*66;
        
        ICE_gear_ratio=TRAN.ratio_gear(gear)*TRAN.ratio_diff; %engine gear ratio
        v_w=ice_w_e/ICE_gear_ratio;  % rotational speed
        em_w_e=EM.gearbox_ratio*v_w; % em speed
        
        
        
            for n1=0:m1 % torque sweep
                Treq=n1*54;

                % look up engine peak torque at that speed
                ice_maxT = interp1(ICE.maxT_w,ICE.maxT_T,ice_w_e);

                % look up  EM peak torque at that speed
                em_maxT = interp1(new_em_maxT_w,new_em_maxT_T,em_w_e);

                %search range

 %               T_upper=min(ice_maxT,Treq+em_maxT);
 %              T_lower=max(Treq-em_maxT,0);

                %discretization
  %             ice_trq=linspace(T_lower,T_upper,20);
  %            em_trq = Treq-ice_trq;
                
                em_trq=linspace(rtl,em_maxT,20);
                ice_trq= Treq-em_trq;


                %% ICE
                % look up enigne fuel rate
                m_equ_ice=interp2(ICE.map_w,ICE.map_T,ICE.map_m_dot',ice_w_e,ice_trq);


                    %% EM
                    % look up table motor efficiency map
                    for n=1:size(em_trq,2) % speed sweep
                    s_chg=interp2(new_em_map_w,new_em_map_T,new_em_map_eff',em_w_e,em_trq(n));
                    s_dischg=s_chg;

                    gamma = (1+sign(em_trq(n)))/2;
                    m_equ_em_no_penalty(n) =(gamma/s_chg + (1-gamma)/s_dischg).*(em_trq(n)*em_w_e)/Qhlv;
                    
                   
                    m_equ_em = f_penalty * m_equ_em_no_penalty;
                    
                    end

                    m_eqv = m_equ_ice + m_equ_em';

%                     plot(m_eqv,'r');
%                     hold on
%                     plot(m_equ_ice,'b');
%                     plot(m_equ_em,'k');

                    [C,I] = min(m_eqv);
                    em_trq_e(n2+1,n1+1,n3)=em_trq(I);
                    
            end   
    end
 end
 
  if gear==1
      em_trq_e_g1=em_trq_e;
    else if gear==2
          em_trq_e_g2=em_trq_e;
          else if gear==3
                  em_trq_e_g3=em_trq_e;
                  else if gear==4
                          em_trq_e_g4=em_trq_e;
                      end
              end
          end
  end
  
          

 
 end

 Treq=linspace(0,20*m1,m1+1);
 ice_w_e=linspace(0,50*m2,m2+1);  
 soc_e=linspace(-1,1,m3);
 
 em_trq_e=em_trq_e_g1;
 figure(1)
 surf(ice_w_e,Treq,em_trq_e(:,:,1)')
 xlabel('engine speed rad/s')
 ylabel('torque request Nm')
 title('EM torque x_soc -1')
 
 figure(2)
 surf(ice_w_e,Treq,em_trq_e(:,:,2)')
 xlabel('engine speed rad/s')
 ylabel('torque request Nm')
 title('EM torque at x_soc -0.5')
 
 figure(3)
 surf(ice_w_e,Treq,em_trq_e(:,:,3)')
 xlabel('engine speed rad/s')
 ylabel('torque request Nm')
 title('EM torque at x_soc 0')
 
 figure(4)
 surf(ice_w_e,Treq,em_trq_e(:,:,4)')
 xlabel('engine speed rad/s')
 ylabel('torque request Nm')
 title('EM torque at x_soc 0.5')
 
 figure(5)
 surf(ice_w_e,Treq,em_trq_e(:,:,5)')
 xlabel('engine speed rad/s')
 ylabel('torque request Nm')
 title('EM torque at x_soc 1')
 
% %  save('ecms_table_n5.mat','Treq','ice_w_e','soc_e','em_trq_e_g1','em_trq_e_g2','em_trq_e_g3','em_trq_e_g4')
 