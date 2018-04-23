function [Accel_part,Force] = Force_Computation(Charge,ElectricField,Length_dom,Pos_cells,Pos_part,Numb_cells,Mass)

    Sum = 0;
    DeltaX = Length_dom/Numb_cells;
    a = zeros(1,length(Pos_cells));
    
            for k = 1:length(Pos_cells)
                a(k) = 1-(abs(Pos_part-Pos_cells(k))/DeltaX);
                    if a(k) < 0
                        a(k) = 0;
                    end
                Sum = Sum+ElectricField(k)*a(1,k);
            end
    
    Force = Charge*Sum*DeltaX;
    Accel_part = Force/Mass;
    
end
                
        
       
    