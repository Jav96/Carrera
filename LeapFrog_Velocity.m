function Vel_part=LeapFrog_Velocity(Vel_part,Deltat,Accel_part,Numb_part,Numb_cells)
         
    %Update Velocity
        for i = 1:Numb_part
            for j = 1:Numb_cells
                Vel_part(i,j)=Vel_part(i,j)+0.5*Deltat*(Accel_part(i,j,1)+Accel_part(i,j,2));
            end
        end
end