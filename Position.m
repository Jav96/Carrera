function Pos_part= Position(Length_dom,Numb_cells,Pos_part)

    Length_cell=Length_dom/Numb_cells;
    [Rows,Columns]=size(Pos_part);
    
    for i=1:Rows
        for j=1:Columns
            dx=Length_cell*rand;
            Pos_part(i,j)=dx+(i-1)*Length_cell;
        end
    end
end
   