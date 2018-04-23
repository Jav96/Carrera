function Index_part = Index_cells(Pos_part,Length_dom,Numb_cells)
    
    Length_cell = Length_dom/Numb_cells;
    
    for i = 1:Numb_cells
        [a,b] = find(Pos_part<=Length_cell*i & Pos_part>=Length_cell*(i-1));
        for j = 1:length(a)
                Index_part(a(j),b(j))=i;
        end
    end
end