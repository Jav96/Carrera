clc;clear
V_mean=2;
Sigma=4;

V_min=V_mean-4*Sigma;
V_max=V_mean+4*Sigma;
Vel_part=(V_max-V_min)*rand+V_min;
if Vel_part<=V_max && Vel_part>=V_min
    fprintf('The velocity is inside the range [V_min;V_max]\n')
    fprintf('Vmin=%f Vel_part=%f V_max=%f\n',V_min,Vel_part,V_max)
else
    fprintf('The velocity does not belong to the range [V_min;V_max]\n')
end