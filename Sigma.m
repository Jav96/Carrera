function Weight=Sigma(Pos_part,Numb_cells,Pos_cells,Numb_part,Index_part,Length_dom)

    DeltaX=Length_dom/Numb_cells;

    for i=1:Numb_cells
        for j=1:Numb_part
            Weight(i,j,1)=1-abs(Pos_part(i,j)-Pos_cells(1,Index_part(i,j)))/DeltaX;
            Weight(i,j,2)=1-abs(Pos_part(i,j)-Pos_cells(1,Index_part(i,j)+1))/DeltaX;
        end
    end
end
                