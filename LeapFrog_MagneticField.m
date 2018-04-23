function [Vel_par,Pos_par,Part_1] = LeapFrog_MagneticField(Pos_part,Vel_part,Mass,Length_dom,Charge,Pos_cells,Numb_part,Numb_cells,Index_part,Deltat,B)

    [~,~,Node_chargeDensity]=Node(Mass,Pos_part,Pos_cells,Numb_cells,Vel_part,Length_dom,Index_part,Charge);
    [~,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);


    for i = 1:Numb_cells
        for j = 1:Numb_part
            
            [Accel_part(i,j),~] = Force_Computation(Charge(i,j),ElectricField,Length_dom,Pos_cells,Pos_part(i,j),Numb_cells,Mass(i,j));
          
            % z-axis
            Vel_par(i,j,1) = Accel_part(i,j)*Deltat+Vel_part(i,j,1);
            Pos_par(i,j,1) = Pos_part(i,j,1)+Vel_par(i,j,1)*Deltat;
            
            if Pos_par(i,j,1) >= Length_dom || Pos_par(i,j,1) < 0
                 x = floor(Pos_par(i,j,1)/Length_dom);
                 Pos_par(i,j,1) = Pos_par(i,j,1)-x*Length_dom;
            end
            
           % x and y axis
           f = @(x)root_vel(x,Vel_part(i,j,2),Vel_part(i,j,3),Deltat,Charge(i,j),Mass(i,j),B);
%             f(1) = @(x)(x(1)-Vel_part(i,j,2))/Deltat-Charge(i,j)*(x(2)+Vel_part(i,j,3))*B/(2*Mass(i,j));
%             f(2) = @(x)(x(2)-Vel_part(i,j,3))/Deltat+Charge(i,j)*(x(1)+Vel_part(i,j,2))*B/(2*Mass(i,j));
            x0 = [2e5,2e5];
            x = fsolve(f,x0);
            
            Vel_par(i,j,2) = x(1);
            Vel_par(i,j,3) = x(2);
            
            Pos_par(i,j,2) = Pos_part(i,j,2)+Vel_par(i,j,2)*Deltat;
            Pos_par(i,j,3) = Pos_part(i,j,3)+Vel_par(i,j,3)*Deltat;

%             f = @(x)(x-Vel_part(i,j,3))/Deltat+Charge(i,j)*(x+Vel_part(i,j,2))*B/(2*Mass(i,j));
%             Vel_par(i,j,3) = fzero(f,2e5);
%             Pos_par(i,j,3) = Pos_part(i,j,1)+Vel_par(i,j,3)*Deltat;
            

            
        end
    end
            Part_1(1) = Pos_par(1,1,1);
            Part_1(2) = Pos_par(1,1,2);
            Part_1(3) = Pos_par(1,1,3);
end
            
            

    