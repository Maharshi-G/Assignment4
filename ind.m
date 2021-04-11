function ind(n1, n2, val)
global G C b   %define global variables
d = size(G,1); % current size of the MNA
xr = d+1;      % new (extera)row/column

% Using an index bigger than the current size,  Matlab automatically 
... increases the size of the matrix:

G(xr,xr) = 0; % add new row/column
C(xr,xr) = 0;
b(xr) = 0;    % add new row

if (n1 ~= 0)
    G(n1,xr) = 1;
    G(xr,n1) = 1;
end
if (n2 ~= 0)
    G(n2,xr) = -1;
    G(xr,n2) = -1;
end

    C(xr,xr) = - val;

end