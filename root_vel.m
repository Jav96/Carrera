function f = root_vel(x,Vel_part1,Vel_part2,Deltat,Charge,Mass,B)
        f(1) = (x(1)-Vel_part1)/Deltat-Charge*(x(2)+Vel_part2)*B/(2*Mass);
        f(2) = (x(2)-Vel_part2)/Deltat+Charge*(x(1)+Vel_part1)*B/(2*Mass);
end
        
       
