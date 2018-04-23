function Figure = Plot(Electric_Field,Node_chargeDensity,Length_dom,Numb_cells,Figure)

    DeltaX = Length_dom/Numb_cells;
    for i = Numb_cells+1
        figure(Figure)
        plot(DeltaX*(i-1),Electric_Field(i))
        hold on
        figure(Figure+1)
        plot(DeltaX*(i-1),Node_chargeDensity(i))
        hold on
    end
    Figure = Figure + 2;
end
