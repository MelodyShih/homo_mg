%
% Antti Hannukainen 4.3.2018 / Otaniemi
%
%----------------------------------------------------------------- 
%
% Assembly routine for triangular P1-elements. Standard "hat" basis
% functions are used. The function computes matrix K and vector F. No 
% boundary conditions are imposed in K or F. 
%
% The matrix K corresponds to finte element discretization of the 
% the bilinear form
%
% \int_\Omega \nabla u^T*A*\nabla v \dx + c*u*v.
%
% where A is piecewise constant over given triangulation such that
% that A_|K = At(K). The vector F corresponds to finte element 
% discretization of the linear functional
%
% \int_\Omega fun*v.
%
% where fun is the load function.
%
% CALLING SYNTAX IS
%
% function [K,F] = assembly_P1(mesh,At,c,fun)
%
%   mesh    = mesh structure representing the triangulation, see inittri.m
%   
%   At      = Nt-vector of the value of coefficient A over mesh. 
%             A_|K = At(K). Scalar At is appropriately extended vector of
%             correct size.
%              
%   c       = scalar coefficient for the L^2 - term in the bilinear form.
%   
%   fun     = function handle fun(x,y). The function has to accept matrix
%             valued input x,y and return matrix value value of the same
%             size such that val = fun(x,y) corresponds to 
%             val(i,j) = fun(x(i,j), y(i,j)). 
%
%     K     = Np x Np sparse matrix corresponding to P1-FE discretization
%             of the bilinear form \int_\Omega \nabla u^T*A*\nabla v \dx +
%             c*u*v. Here Np is the number of nodes in mesh. 
%
%     F     = Np x 1 full vector corresponding to P1-FE discretization
%             of the load functional \int_\Omega fun*v.
%
% SEE : example_poisson.m 
%

function [K,F] = assembly_P1(mesh,At,c,fun)

Ndof = size(mesh.p,2);

% Compute vectorised affine mapping.
map = affine_tri(mesh);

X = [1/2 1/2 0 ; 0 1/2 0];
W = [1/6 1/6 1/6]';

gX{1} = bsxfun(@plus,map.Ax*X,map.bx);
gX{2} = bsxfun(@plus,map.Ay*X,map.by);

% Evaluate the load
fval = fun(gX{1},gX{2});

% Evaluate the coefficient
if( isscalar(At))
   Aval = repmat(At,size(gX{1}));
else
   Aval = repmat(At(:),1,3); 
end

% Define reference basisfunctions
% and their derivatives
L{1} = 1-X(1,:)-X(2,:);
L{2} = X(1,:);
L{3} = X(2,:);

Nip = size(X,2);
Nt = size(mesh.t,2);

dL{1} = [ -ones(1,Nip); -ones(1,Nip) ];
dL{2} = [  ones(1,Nip); zeros(1,Nip) ];
dL{3} = [ zeros(1,Nip);  ones(1,Nip) ];

% Preallocate for speed
ffind = [];
ff = zeros(3,Nt);

iind = [];
jind = [];
kk = zeros(9,Nt);

mind = 1;

% The assembly loop
for i=1:3
    
    Li = repmat(L{i},Nt,1);
    
    dLix = map.Px*dL{i};
    dLiy = map.Py*dL{i};
    
    ff(i,:) = (fval.*Li)*W.*abs(map.detA);
    
    for j=1:3
        
        Lj = repmat(L{j},Nt,1);
        
        
        dLjx = map.Px*dL{j};
        dLjy = map.Py*dL{j};
               
        % Keep track of indeces 
        iind = [iind i] ; jind = [jind j];
        
        kk(mind,:) = ( Aval.*(dLix.*dLjx + dLiy.*dLjy) + c*(Li.*Lj) ) *W.*abs(map.detA);
        mind = mind+1;
    end
end

K = sparse(mesh.t(iind,:),mesh.t(jind,:),kk,Ndof,Ndof);
F = full(sparse(mesh.t,ones(size(mesh.t)),ff,Ndof,1));

