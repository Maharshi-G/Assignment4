function ccvs(nd1,nd2,ni1,ni2,val)
% Adds the stamp of a current controlled voltage source 
% G matrix in circuit representation.
%   ni1 -------o            |----------o nd1
%              |            |
%              |           /+\
%          Ii| |          /   \    Vni1 - Vni2 = val*(Ii)
%            | |          \   /
%           \ /|           \-/ 
%              |            |
%   ni2 -------o            |----------o nd2
global G C b;
d = size(G,1);
b(d+1) = 0; 
b(d+2) = 0; 
G( d+1, d+1) = 0; 
G(d+2,d+2) = 0;
C( d+1, d+1) = 0;
C(d+2,d+2) = 0;
if (ni1 ~= 0)
    G( d+1,ni1) = 1;
    G(ni1, d+1) = 1;
end
if (ni2 ~= 0)
    G( d+1,ni2) = -1;
    G(ni2, d+1) = -1;
end
if (nd1 ~= 0)
    G(d+2,nd1) = 1;
    G(nd1,d+2) = 1;
end
if (nd2 ~= 0)
    G(d+2,nd2) = -1;
    G(nd2,d+2) = -1;
end
G(d+2, d+1) = -val;
%END