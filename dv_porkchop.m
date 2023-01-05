function Z = dv_porkchop(X, Y, p1, p2, grid,mu)
% dv_porkchop: Calculates a grid of delta-v values
%
% INPUTS:
% X, Y = Grids of departure and arrival dates
% p1, p2 = Departure and arrival planets
% grid = Function handle to compute delta-v value
% mu = Gravitational parameter of the central body
%
% OUTPUTS:
% Z = Matrix of delta-v values
%
% USAGE:
% Z = dv_porkchop(X, Y, p1, p2, grid, mu)
% Authors
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it

i = 0;
% Loop through X and Y and compute delta-v values
while(i<size(X,1))
    i = i+1;
    j = 0;
    while(j<size(X,2))
        j = j+1;
        Z(i,j) = grid(X(i,j), Y(i,j), p1, p2, mu);
    end
end

end