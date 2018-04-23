function [Pos_par,Accel_part]=LeapFrog_Position(Pos_part,Vel_part,Mass,Length_dom,Charge,Pos_cells,Numb_part,Deltat,Numb_cells,Index_part)
     
    [~,~,Node_chargeDensity]=Node(Mass,Pos_part,Pos_cells,Numb_cells,Vel_part,Length_dom,Index_part,Charge);
    [~,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);
    [Accel_part,~] = Force_Computation(Charge,ElectricField,Length_dom,Pos_cells,Pos_part,Numb_cells,Numb_part,Mass);
    
    %Update Position
            Pos_par=Pos_part+Vel_part*Deltat+0.5*Accel_part*Deltat^2;
                if Pos_par >= Length_dom
                    x = floor(Pos_par/Length_dom);
                    Pos_par = Pos_par-x*Length_dom;
                end

end