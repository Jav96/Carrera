function [Pos_par,Vel_par] = LeapFrog(Pos_part,Vel_part,Mass,Length_dom,Charge,Pos_cells,Numb_part,Numb_cells,Deltat,u)
    
        Node_chargeDensity=ChargeDensity(Pos_cells,Pos_part,Charge,Numb_part,Numb_cells,Length_dom);
        [~,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);
    
        if u == 1   
            for i = 1:Numb_cells
                for j = 1:Numb_part
                    [Accel_part(i,j),~] = Force_Computation(Charge(i,j),ElectricField,Length_dom,Pos_cells,Pos_part(i,j),Numb_cells,Mass(i,j));
                    Vel_par(i,j) = Vel_part(i,j)-Accel_part(i,j)*Deltat/2;
                    Vel_par(i,j) = Accel_part(i,j)*Deltat+Vel_par(i,j);
                    Pos_par(i,j)=Vel_par(i,j)*Deltat+Pos_part(i,j);
        %                 if Pos_par(i,j) >= Length_dom || Pos_par(i,j) < 0
        %                     x = floor(Pos_par(i,j)/Length_dom);
        %                     Pos_par(i,j) = Pos_par(i,j)-x*Length_dom;
        %                 end
                end
            end
        else
            for i = 1:Numb_cells
                for j = 1:Numb_part
                    [Accel_part(i,j),~] = Force_Computation(Charge(i,j),ElectricField,Length_dom,Pos_cells,Pos_part(i,j),Numb_cells,Mass(i,j));
                    Vel_par(i,j) = Accel_part(i,j)*Deltat+Vel_part(i,j);
                    Pos_par(i,j)=Vel_par(i,j)*Deltat+Pos_part(i,j);
        %                 if Pos_par(i,j) >= Length_dom || Pos_par(i,j) < 0
        %                     x = floor(Pos_par(i,j)/Length_dom);
        %                     Pos_par(i,j) = Pos_par(i,j)-x*Length_dom;
        %                 end
                end
        end
    end

    




