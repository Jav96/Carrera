%% Slabs of electrons displaced < DeltaX/2 from static slabs of ions (Check PIC code)
clear all;
clc;
 %% Variables
 
    Length_dom = 100;
    
 %% Initialized Data
 
    Numb_cells=input('Select the number of desired cells:');
    while mod(Numb_cells,1)~=0
        Numb_cells=input('Select an integer value please:');
    end
    
    Numb_part = 2;
    
    DeltaX = Length_dom/Numb_cells;
    DeltaX_electrons = rand*DeltaX/2;
    Time = input('Select time of computation:');
    Steps = input('Select the number of desired steps:');
    Deltat = Time/Steps;
    for i = 1:Steps
        Vector_time(i) = Deltat*i;
    end

    for i = 1:Numb_cells+1
        Pos_cells(i) = Length_dom*(i-1)/Numb_cells;
    end
    
    Mass = zeros(Numb_cells,Numb_part);
    Charge = zeros(Numb_cells,Numb_part);
    Pos_part = zeros(Numb_cells,Numb_part);
    Vel_part = zeros(Numb_cells,Numb_part);
    

%     for i = 1:Numb_cells
%         for j = 1:Numb_part
%         a=rand;
%             if a<=0.5
%                 Mass(i,j) = 9.1e-31;
%                 Charge(i,j) = -1.6e-19;
%                 Pos_part(i,j) = Pos_cells(i)+DeltaX_electrons;
%                 Vel_part(i,j) = 0;
%             else
%                 Mass(i,j) = 2.18; %mi >>> me
%                 Charge(i,j) = 1.6e-19;
%                 Pos_part(i,j) = Pos_cells(i);
%                 Vel_part(i,j) = 0;
%             end
%         end
%     end      

%Set Numb_cells = 2
for j = 1:Numb_cells
for i = 1:Numb_part
    if mod(j,2) == 0
        Mass(j,i) = 9.1e-31;
        Charge(j,i) = -1.6e-19;
        Pos_part(j,i) = Pos_cells(j)+DeltaX_electrons;
    else
        Pos_part(j,i) = Pos_cells(j);  
        Charge(j,i) = 1.6e-19;
        Mass(j,i) = 2.18;
    end
end
end

    for h = 1:Steps
        [Pos_part,Vel_part] = LeapFrog(Pos_part,Vel_part,Mass,Length_dom,Charge,Pos_cells,Numb_part,Numb_cells,Deltat,h);
        for i = 1:Numb_part
            if Mass(2,i) == 9.1e-31
                 Pos_par(h) = Pos_part(2,i);
                 Vel_par(h) = Vel_part(2,i);
                 break
            end
        end   
    end
    
    figure(1)
    plot(Vector_time,Pos_par);
    title('Electron Position vs Time')
    figure(2)
    plot(Vector_time,Vel_par);
    title('Electron Velocity vs Time')
    Figure = 3;
    
    Node_chargeDensity = ChargeDensity(Pos_cells,Pos_part,Charge,Numb_part,Numb_cells,Length_dom);
    [Phi,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);
    
        figure(Figure)
        plot(Pos_cells,ElectricField)
        title('ElectricField')
        figure(Figure+1)
        plot(Pos_cells,Node_chargeDensity)
        title('Node Charge Density')
        figure(Figure+2)
        plot(Pos_cells,Phi)
        title('Electric Potential')
        Figure = Figure + 3;
        figure(Figure)
        plot3(Vector_time,Vel_par,Pos_par);
        grid on
        xlabel('Time')
        ylabel('Velocity')
        zlabel('Position')
    
    
%     figure(Figure)
%     plot(Pos_part(1,1)+1,0.075,'.','LineWidth',2);
    
