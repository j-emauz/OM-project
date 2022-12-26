function [Z,X,Y] = dv_totalcost(X1, Y1_GA,X2_GA,Y2, p1, p2,NEO,mu)
% Z = dv_calcgrid(X, Y, p1, p2, mu)
%   Detailed explanation goes here



i = 0;
while(i<size(X1,1))
    i = i+1;
    j = 0;
    while(j<size(X1,1))
        j = j+1;
        % Z(i,j) = grid(X(i,j), Y(i,j), p1, p2, mu);
        Z(i,j) = dv_arc1(X1(i,j), Y1_GA(i,j), p1, p2, mu);
%         Z(i,j) = total_cost(X1(i,j), Y1_GA(i,j),X2_GA(m,t),Y2(m,t), p1, p2,NEO, mu);
        % was Z(i,j) = total_cost(X(i,j), Y(i,j)+ Y_GA(i,j),X(i,j)+X_GA(i,j),Y(i,j), p1, p2,NEO, mu);
        % function [dv_tot,VF,v_2, car_2,t1,dt] = total_cost(t1, t2,t3,t4, p1, p2,NEO, mu)

        % Saturn    
        % [X, Y] = meshgrid(tspan_dept, tspan_GA);
        % Z = dv_porkchop(X, Y, p1, p2, @dv_arc1,mu_S);

        % Asteroid
        % [X, Y] = meshgrid(tspan_GA, tspan_arrt);
        % Z = dv_porkchop(X, Y, p2,86, @dv_arcNEO,mu_S);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pt 2
m = 0;
while(m<size(X2_GA,1))
    m = m+1;
    t = 0;
    while(t<size(X2_GA,1))
        t = t+1;
        % Z(i,j) = grid(X(i,j), Y(i,j), p1, p2, mu);
        Z(i+m,j+t)= Z(i,j)+dv_arcNEO(X2_GA(i,j), Y2(i,j), p2, NEO, mu);
%         Z(i,j) = total_cost(X1(i,j), Y1_GA(i,j),X2_GA(m,t),Y2(m,t), p1, p2,NEO, mu);
        % was Z(i,j) = total_cost(X(i,j), Y(i,j)+ Y_GA(i,j),X(i,j)+X_GA(i,j),Y(i,j), p1, p2,NEO, mu);
        % function [dv_tot,VF,v_2, car_2,t1,dt] = total_cost(t1, t2,t3,t4, p1, p2,NEO, mu)
        X(i+m,j+t)=X1(i,j)+X2_GA(m,t);
        Y(i+m,j+t)=Y1_GA(i,j)+Y2(m,t);
        % Saturn    
        % [X, Y] = meshgrid(tspan_dept, tspan_GA);
        % Z = dv_porkchop(X, Y, p1, p2, @dv_arc1,mu_S);

        % Asteroid
        % [X, Y] = meshgrid(tspan_GA, tspan_arrt);
        % Z = dv_porkchop(X, Y, p2,86, @dv_arcNEO,mu_S);
    end
end


end