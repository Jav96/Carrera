function [Phi,ElectricField]=Poisson(Numb_cells,Length_dom,Node_chargeDensity,Pos_cells)

    Epsilon=8.8541878176e-12;
    DeltaX=Length_dom/Numb_cells;
    Case = 1; %1 == Dirichlet condition // 2 == Robin condition (Left Neumann, Right Dirichlet)
              %3 == Robin Conditions (Right Neumann, Left Dirichlet)
    
    if Case == 1

    MatrixA=zeros(length(Pos_cells)-2);
    
    % Boundary Conditions (Alfa==>Beginning of domain;Beta==>End of domain)
    Alpha = 0;
    Beta = 0;
    
        if Numb_cells == 1
            Phi = [Alpha,Beta];
        else
        
    for i=1:Numb_cells-2
            MatrixA(i,i) = -2/DeltaX^2;
            MatrixA(i+1,i) = 1/DeltaX^2;
            MatrixA(i,i+1) = 1/DeltaX^2;
    end
        MatrixA(Numb_cells-1,Numb_cells-1) = -2/DeltaX^2;


%     [Size,~]=size(MatrixA);
%     MatrixA=MatrixA(1:Size-1,1:Size-1);
    ColumnP=Node_chargeDensity(1,2:length(Node_chargeDensity)-1)/Epsilon;
    ColumnP(1)=ColumnP(1)-Alpha/DeltaX^2;
    ColumnP(length(ColumnP))=ColumnP(length(ColumnP))-Beta/DeltaX^2;
    Phi=MatrixA\ColumnP';
    Phi=[Alpha,Phi',Beta];
        end
    
    ElectricField=zeros(1,length(Pos_cells));
    for i=1:length(Pos_cells)
        if i == 1 
            ElectricField(i)=Phi(i+1)-Phi(i)/DeltaX;
        elseif i == length(Pos_cells)
            ElectricField(i)=Phi(i)-Phi(i-1)/DeltaX;
        else
            ElectricField(i)=Phi(i+1)-Phi(i-1)/2*DeltaX;
        end
    end
    elseif Case == 2 
        
        % Boundary Conditions (Alpha==>Beginning of domain(Neumann);Beta==>End of domain(Dirichlet)
        Alpha = 0;
        Beta = 0;
    
        MatrixA = zeros(Numb_cells+1);
        ColumnP = zeros(Numb_cells+1,1);
        
        MatrixA(1,1)=-1/DeltaX^2;
        MatrixA(1,2)=2/DeltaX^2;
        %ColumnP(2)=0.5*(Node_chargeDensity(1)/Epsilon+2*DeltaX*Alpha/DeltaX^2);
        for i = 2:Numb_cells
            MatrixA(i,i-1) = 1/DeltaX^2;
            MatrixA(i,i) = -2/DeltaX^2;
            MatrixA(i,i+1) = 1/DeltaX^2;
        end
            MatrixA(Numb_cells+1,Numb_cells) = 1/DeltaX^2;
            MatrixA(Numb_cells+1,Numb_cells+1) = -2/DeltaX^2;
            for i = 1:Numb_cells+1
                ColumnP(i,1) = -Node_chargeDensity(i)/Epsilon;
            end
                ColumnP(1,1) = -Node_chargeDensity(1)/(Epsilon*2)+Alpha/DeltaX;
                ColumnP(Numb_cells+1,1) = -Node_chargeDensity(Numb_cells+1)/(Epsilon)-Beta/DeltaX^2;
%         MatrixA(Numb_cells+3,Numb_cells+2)=P;
%         MatrixA(Numb_cells+3,Numb_cells+3)=-P;
%         MatrixA(1,2)=P;
%         ColumnP(Numb_cells+2)=0.5*(Node_chargeDensity(Numb_cells+1)/Epsilon+2*DeltaX*Beta/DeltaX^2);
%         ColumnP(1)=2*Alpha/DeltaX;
%         ColumnP(Numb_cells+3)=2*Beta/DeltaX;
        
        Phi=MatrixA\(ColumnP);
%         Phi=Phi(2:length(Phi)-1);
        ElectricField=zeros(1,length(Pos_cells));
    for i=1:length(Pos_cells)
        if i == 1 
            ElectricField(i)=Phi(i+1)-Phi(i)/DeltaX;
        elseif i == length(Pos_cells)
            ElectricField(i)=Phi(i)-Phi(i-1)/DeltaX;
        else
            ElectricField(i)=Phi(i+1)-Phi(i-1)/DeltaX;
        end
    end
    end
        
        
        
            
   end 
    
    