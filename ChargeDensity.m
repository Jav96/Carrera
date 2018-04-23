function Node_chargeDensity = ChargeDensity (Pos_cells,Pos_part,Charge,Numb_part,Numb_cells,Length_dom)

    DeltaX = Length_dom/Numb_cells; 
    
    for i = 1:length(Pos_cells)
        Node_chargeDensity(i) = 0;
        for j = 1:Numb_cells
            for k = 1:Numb_part
                a = (1-abs(Pos_cells(i)-Pos_part(j,1))/DeltaX);
                if a < 0
                    a = 0;
                else
                    Node_chargeDensity(i) = Node_chargeDensity(i) + a*Charge(j,1);
                end
            end
        end
    end
end