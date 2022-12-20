function Z = dv_calcgrid(X, Y, p1, p2, mu)
% Z = dv_calcgrid(X, Y, p1, p2, mu)
%   Detailed explanation goes here

i = 0;
while(i<size(X,1))
    i = i+1;
    j = 0;
    while(j<size(X,2))
        j = j+1;
        Z(i,j) = dv_calc(X(i,j), Y(i,j), p1, p2, mu);
    end
end

end