function [Vel_part,f_star] = Velocity (T,Mass,V_mean,k,Numb_part)
 
    Sigma=sqrt(T*k/Mass);
    V_min=V_mean-4*Sigma;
    V_max=V_mean+4*Sigma;
        
    Vel_part=(V_max-V_min)*rand+V_min;
    f_star=exp((-Mass*((Vel_part-V_mean)^2))/(2*k*T));
    f_prime=rand;
        
        while f_star>f_prime
             Vel_part=(V_max-V_min)*rand+V_min;
             f_star=exp(-Mass*(Vel_part-V_mean)^2/(2*k*T));
             f_prime=rand;
        end
        
end