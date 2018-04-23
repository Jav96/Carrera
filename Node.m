function [Node_density,Node_velocity,Node_chargeDensity,Node_currentDensity]=Node(Mass,Pos_part,Pos_cells,Numb_cells,Vel_part,Length_dom,Index_part,Charge,Numb_part)
    
    DeltaV=1;
    DeltaX=Length_dom/Numb_cells;
    
    if Numb_cells == 1
        for i = 1:2
            Sum=zeros(1,length(Pos_cells));
            Sum2=zeros(1,length(Pos_cells));
            Sum3=zeros(1,length(Pos_cells));
            Sum4=zeros(1,length(Pos_cells));
            for j = 1:Numb_part
                Sum(i) = Sum(i)+Mass(1,j)*(1-abs(Pos_cells(i)-Pos_part(1,j)))/DeltaX;
                Sum2(i) = Sum2(i)+Mass(1,j)*Vel_part(1,j)*(1-abs(Pos_cells(i)-Pos_part(1,j)))/DeltaX;
                Sum3(i) = Sum3(i)+Charge(1,j)*(1-abs(Pos_cells(i)-Pos_part(1,j)))/DeltaX;
                Sum4(i) = Sum4(i)+Charge(1,j)*Vel_part(1,j)*(1-abs(Pos_cells(i)-Pos_part(1,j)))/DeltaX;
            end
            Node_density(i) = Sum(i)/DeltaV;
            Node_velocity(i) = Sum2(i)/(DeltaV*Node_density(i));
            Node_chargeDensity(i) = Sum3(i);
            Node_currentDensity(i) = Sum4(i);
        end
        
    else        
    
    for i = 1:length(Pos_cells)
        if i==1
            [a,b]=find(Index_part==1);
            [d,e]=find(Index_part==i);
        elseif i==length(Pos_cells)
            [a,b]=find(Index_part==i-1);
            [d,e]=find(Index_part==i);
        else
            [a,b]=find(Index_part==i-1); 
            [d,e]=find(Index_part==i);
        end
        c=[a b];
        f=[d e];
        
        Sum=zeros(1,length(Pos_cells));
        Sum2=zeros(1,length(Pos_cells));
        Sum3=zeros(1,length(Pos_cells));
        Sum4=zeros(1,length(Pos_cells));
        for j=1:length(a)
            if i==1
                Sum(i)=Sum(i)+Mass(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum2(i)=Sum2(i)+Mass(c(j,1),c(j,2))*Vel_part(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum3(i)=Sum3(i)+Charge(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum4(i)=Sum4(i)+Charge(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX))*Vel_part(c(j,1),c(j,2));
            elseif i==length(Pos_cells)
                Sum(i)=Sum(i)+Mass(c(j,1),c(j,2))*(1-abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX);
                Sum2(i)=Sum2(i)+Mass(c(j,1),c(j,2))*Vel_part(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum3(i)=Sum3(i)+Charge(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum4(i)=Sum4(i)+Charge(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX))*Vel_part(c(j,1),c(j,2)); 
            else
                Sum(i)=Sum(i)+Mass(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum2(i)=Sum2(i)+Mass(c(j,1),c(j,2))*Vel_part(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum3(i)=Sum3(i)+Charge(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX));
                Sum4(i)=Sum4(i)+Charge(c(j,1),c(j,2))*(1-(abs(Pos_part(c(j,1),c(j,2))-Pos_cells(1,i))/DeltaX))*Vel_part(c(j,1),c(j,2));    
            end
        end
        for j=1:length(d)
            if i~=1 && i~=length(d)
                Sum(i)=Sum(i)+Mass(f(j,1),f(j,2))*(1-(abs(Pos_part(f(j,1),f(j,2))-Pos_cells(1,i)/DeltaX)));
                Sum2(i)=Sum2(i)+Mass(f(j,1),f(j,2))*Vel_part(f(j,1),f(j,2))*(1-(abs(Pos_part(f(j,1),f(j,2))-Pos_cells(1,i))/DeltaX));
                Sum3(i)=Sum3(i)+Charge(f(j,1),f(j,2))*(1-(abs(Pos_part(f(j,1),f(j,2))-Pos_cells(1,i))/DeltaX));
                Sum4(i)=Sum4(i)+Charge(f(j,1),f(j,2))*(1-(abs(Pos_part(f(j,1),f(j,2))-Pos_cells(1,i))/DeltaX))*Vel_part(f(j,1),f(j,2));
            end
        end
        Node_density(i) = Sum(i)/DeltaV;
        Node_velocity(i) = Sum2(i)/(DeltaV*Node_density(i));
        Node_chargeDensity(i) = Sum3(i);
        Node_currentDensity(i) = Sum4(i);
    end
    end
end