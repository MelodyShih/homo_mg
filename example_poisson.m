%
% Antti Hannukainen 4.3.2018 / Otaniemi
%
%----------------------------------------------------------------- 
%
% Solve the Poisson equation on domain (0,r)x(0,r) with zero dirichlet
% boundary condition and load function 1.  Finite element method with
% triangular P1-elements is used.
%
% requires util-folder in the path
%


% make the mesh (use Nref refinement steps. Domain is (0,r)x(0,r) ).

r = 2;
Nref = 7; 
mesh = make_Ur_mesh(r,Nref);

% assembly the FE-matrix and load

[K,b] = assembly_P1(mesh,1,0, @(x,y)(ones(size(x))) );

% set zero dirichlet BC.
in = mesh.in;
x = zeros( size(mesh.p,2),1);
x(in) = K(in,in)\b(in);

% Plot the solution.
X = mesh.p(1,:);
Y = mesh.p(2,:);
t = mesh.t;

figure;
P = patch(X(t),Y(t),x(t),x(t));
set(P,'edgecolor','none');
