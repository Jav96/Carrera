%% Right Handside Polarized%%
clear all;clc;close all

%% Variables

Length_dom = 5;
Electric_Field_RHS = 1e-8;
Frequency_RHS = 10;
Collisional_Freq = 0;

Numb_cells=input('Select the number of desired cells:');
    while mod(Numb_cells,1)~=0
        Numb_cells=input('Select an integer value please:');
    end
    
Time = input('Select time of computation:');
Steps = input('Select the number of desired steps:');
Deltat = Time/Steps;
DeltaX = Length_dom/Numb_cells;


for i=1:Numb_cells+1
    Pos_cells(i)=Length_dom*(i-1)/Numb_cells;
end

%Magnetic Field z will follow this equation Bz = 0.05z^2 + C; dBz/dz = 0.1z
% for i = 1:Numb_cells+1
%     Bz(i) = 0.05*Pos_cells(i)^2;
%     dBz(i) = 0.1*Pos_cells(i);
% end

Mass = 9.1e-31;
Charge = -1.6e-19*ones(length(Pos_cells),1);



%% Solving equations %% 

%Initial conditions%
x0 = zeros(length(Pos_cells),5);
for i = 1:length(Pos_cells)
    x0(i,1) = 0; %Perpendicular Velocity
    x0(i,2) = 2*pi*i/length(Pos_cells); %Electron phase
    x0(i,3) = 0; %Axial Velocity
    x0(i,4) = Pos_cells(i); %Position
    x0(i,5) = 0; %RHS phase;
end

for i = 1:Steps
    Node_chargeDensity = ChargeDensity(Pos_cells,x0(:,4),Charge,1,Numb_cells,Length_dom);
    [~,Electric_Field] = Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);
    Tspan = [0 Deltat];
    %Solve equation for every particle
    for j = 1:length(Pos_cells)
        Bz = 0.05*x0(j,4)^2;
        dBz = 0.1*x0(j,4);
        if x0(j,4) == 0
            Index = 2;
        else
            Index = find(x0(j,4)<=Pos_cells,1,'first');
        end
        E_part = Electric_Field(Index)*(x0(j,4)-Pos_cells(Index-1))/DeltaX+Electric_Field(Index-1)*(1-((x0(j,4)-Pos_cells(Index-1))/DeltaX));
        [t,x] = ode45(@(t,x) Equation(x,Mass,Charge(j),Electric_Field_RHS,Frequency_RHS,Collisional_Freq,Bz,dBz,Length_dom,E_part),Tspan,x0(j,:)); 
    %Collect data
    if i == 1
        part(j).Time = t;
        part(j).Vel_perp = x(:,1);
        part(j).Beta = x(:,2);
        part(j).Vel_axi = x(:,3);
        part(j).Position = x(:,4);
        part(j).Alpha = x(:,5);
    else
        part(j).Time = [part(j).Time;part(j).Time + t];
        part(j).Vel_perp = [part(j).Vel_perp;x(:,1)];
        part(j).Beta = [part(j).Beta;x(:,2)];
        part(j).Vel_axi = [part(j).Vel_axi;x(:,3)];
        part(j).Position = [part(j).Position;x(:,4)];
        part(j).Alpha = [part(j).Vel_Alpha;x(:,5)];
    end
    %Update initial conditions
    x0(j,:) = x(length(x),:);
    end
end




