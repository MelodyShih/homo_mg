%
% Antti Hannukainen 4.3.2018 / Otaniemi
%
%----------------------------------------------------------------- 
%
% Solve the Poisson equation (A \nabla u,\nabla v) = (f,v) on domain 
% (0,r)x(0,r) with zero dirichlet boundary condition and load function 
% f = 1. The coefficient function A in each cell has value 1 / 9 
% with 50 / 50 probability.
% 
% Finite element method with triangular P1-elements is used.
%
close all;
clear all;

  r = 50;  % size of the domain Ur = (0,r)x(0,r)
Nref = 3;  % number or refinements for the FE - mesh

% generate the mesh
[mesh,t2c] = make_Ur_mesh(r,Nref);

% generate coefficient corresponding to A.
At = make_1_9_cell_At(r, t2c);

% assembly 
[K,b] = assembly_P1(mesh,At,0, @(x,y)(ones(size(x))) );

% set zero dirichlet BC.
bind = mesh.bn;
iind = setdiff(1:size(mesh.p,2),bind);

x = zeros( size(mesh.p,2),1);
x(iind) = K(iind,iind)\b(iind);

% Plot the solution.
X = mesh.p(1,:);
Y = mesh.p(2,:);
t = mesh.t;

figure;
title(['Solution u, r =',num2str(r)]);
P = patch(X(t),Y(t),x(t),x(t));
set(P,'edgecolor','none');

% Plot the coefficient field
figure;
title(['Coefficient field, r =',num2str(r)]);
P = patch(X(t),Y(t),repmat(At(:)',3,1),repmat(At(:)',3,1));
set(P,'edgecolor','none');

