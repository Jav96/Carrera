%% Main %%

clc;clear all;close all

%% Variables

B = 5e-4;
Temp = 10000;
V_mean = 25000;
Length_dom = 5;
k = 1.38e-23;

%% Initialize data

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
    
    %Numb_part=round(10*rand+5);
    Numb_part=20;
    Vel_partz = zeros(Numb_cells,Numb_part);
    Vel_partx = zeros(Numb_cells,Numb_part);
    Vel_party = zeros(Numb_cells,Numb_part);
    
    Pos_partz = zeros(Numb_cells,Numb_part);
    Pos_partx = zeros(Numb_cells,Numb_part);
    Pos_party = zeros(Numb_cells,Numb_part);
    
    Mass = zeros(Numb_cells,Numb_part);
    Charge = zeros(Numb_cells,Numb_part);
    
    
for i = 1:Numb_cells
    for j = 1:Numb_part
        a=rand;
        if a<=1
            Mass(i,j) = 9.1e-31;
            Charge(i,j) = -1.6e-19;
        else
            Mass(i,j) = 2.18e-25;
            Charge(i,j) = 1.6e-19;
        end
        [Vel_partz(i,j),f_star(i,j)] = Velocity(Temp,Mass(i,j),V_mean,k,Numb_part);
        [Vel_partx(i,j),~] = Velocity(Temp,Mass(i,j),V_mean,k,Numb_part);
        [Vel_party(i,j),~] = Velocity(Temp,Mass(i,j),V_mean,k,Numb_part);
        
    end
end

    Pos_partz = Position(Length_dom,Numb_cells,Pos_partz);
    Pos_partx = Position(Length_dom,Numb_cells,Pos_partx);
    Pos_party = Position(Length_dom,Numb_cells,Pos_party);
    
    for i = 1:3
        for j = 1:Numb_part
            for k = 1:Numb_cells
                if i == 1
                    Pos_part(k,j,i) = Pos_partz(k,j);
                    Vel_part(k,j,i) = Vel_partz(k,j);
                elseif i == 2
                    Pos_part(k,j,i) = Pos_partx(k,j);
                    Vel_part(k,j,i) = Pos_partx(k,j);
                else
                    Pos_part(k,j,i) = Pos_party(k,j);
                    Vel_part(k,j,i) = Pos_party(k,j);
                end
            end
        end
    end
    Index_part = Index_cells(Pos_part(:,:,1),Length_dom,Numb_cells);

%% Check distribution
Figure=1;
figure(Figure)
plot(Vel_partz(1,:),f_star(1,:),'.')
Figure=Figure+1;

% for i=1:Numb_cells+1
%     figure(Figure);
%     plot([(Length_dom/Numb_cells)*(i-1)+1,(Length_dom/Numb_cells)*(i-1)+1],[-0.05,0.05],'k','LineWidth',1.1)
%     hold on
% end
% 
% figure(Figure);
% plot([1,Length_dom+1],[0,0],'k','LineWidth',2)
% axis([0 Length_dom+2 -1 1])
% Figure = Figure + 1; 

%  for h = 1:Steps
%          [Vel_part,Pos_part,~] = LeapFrog_MagneticField(Pos_part,Vel_part,Mass,Length_dom,Charge,Pos_cells,Numb_part,Numb_cells,Index_part,Deltat,B);
%          Part_1(h,1) = Pos_part(1,1,1);
%          Part_1(h,2) = Pos_part(1,1,2);
%          Part_1(h,3) = Pos_part(1,1,3);
%  end

for h = 1:Steps
    
    if h==1
    Node_chargeDensity = ChargeDensity(Pos_cells,Pos_partz,Charge,Numb_part,Numb_cells,Length_dom);
    [Phi,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);
    %check
     Sum=sum(Pos_cells.*Node_chargeDensity)*DeltaX;
     Sum2=sum(sum(Charge.*Pos_partz));
     
    
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
    end
    [Pos_partz,Vel_partz] = LeapFrog(Pos_partz,Vel_partz,Mass,Length_dom,Charge,Pos_cells,Numb_part,Numb_cells,Deltat,i);
    
%     figure(Figure)
%     plot(Pos_part(1,1)+1,0.075,'.','LineWidth',2);
    
end
    Node_chargeDensity = ChargeDensity(Pos_cells,Pos_partz,Charge,Numb_part,Numb_cells,Length_dom);
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

% Index_part2 = Index_cells(Pos_part,Length_dom,Numb_cells);
  
 %Weight=Sigma(Pos_part,Numb_cells,Pos_cells,Numb_part,Index_part,Length_dom);
% [Node_density,Node_velocity,Node_chargeDensity]=Node(Mass,Pos_part,Pos_cells,Numb_cells,Vel_part,Length_dom,Index_part,Charge);
% 
% [Phi,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells);
% for i = 1:Numb_cells
%     for j = 1:Numb_part
% [Accel_part(i,j),~] = Force_Computation(Charge(i,j),ElectricField,Length_dom,Pos_cells,Pos_part(i,j),Numb_cells,Mass(i,j));
%     end
% end


% figure(Figure);
% Plot(Numb_cells,Pos_part)
    
 


