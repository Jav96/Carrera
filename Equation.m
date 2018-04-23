function Derivatives = Equation(x,Mass,Charge,Electric_Field_RHS,Frequency_RHS,Collisional_Freq,Bz,dBz,Length_dom,E_part)

%             if x(4) >= Length_dom || x(4) < 0
%         y = floor(x(4)/Length_dom);      %Ensure that if a particle exits, it enters from the beginning
%         x(4) = x(4)-y*Length_dom;
%     end
%     
%     if x(5) >= 2*pi
%         y = floor(x(5)/2*pi);
%         x(5) = x(5)-y*2*pi;
%     end
%                                         %Phases not great that 2pi
%     if x(2) >= 2*pi
%         y = floor(x(2)/2*pi);
%         x(2) = x(2)-y*2*pi;
%     end
% Differential equations with respect to time to solve
        Derivatives = zeros(5,1);
        Derivatives(1) = -Charge*Electric_Field_RHS/Mass*cos(x(5)-x(2))+x(3)*x(1)*dBz/(2*Bz);
        Derivatives(2) = Charge*Bz/Mass - Charge*Electric_Field_RHS*sin(x(5)-x(2))/(Mass*x(1));
        Derivatives(3) = -Charge*E_part/Mass-x(1)^2*dBz/(2*Bz)-Collisional_Freq*x(3);
        Derivatives(4) = x(3);
        Derivatives(5) = Frequency_RHS;
        
        
end