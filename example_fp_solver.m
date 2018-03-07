%
% Antti Hannukainen 4.3.2018 / Otaniemi
%
%----------------------------------------------------------------- 
%
% Solve the Poisson equation (A \nabla u,\nabla v) = (f,v) on domain 
% (0,r)x(0,r) with zero dirichlet boundary condition and load function 
% f = 1. The coefficient function A in each cell has value 1 / 9 
% with 50 / 50 probability. Finite element method with triangular P1-elements is used.
%
% The linear system is solved using the fixed point iteration from Th 1.1.
% The parameter Lambda and the size of the domain r are both fixed. 
% 
% requires util-folder in the path
%

close all;
clear all;

   r = 20;  % size of the domain Ur = (0,r)x(0,r)
Nref = 2;  % number or refinements for the FE - mesh
   L = 0.1; % Lambda parameter for iterative solver.
   
% generate mesh
[mesh,t2c] = make_Ur_mesh(r,Nref);

% generate random pwc. on each cell of Ur.
At = make_1_9_cell_At(r, t2c);

% The corresponding homogenised parameter is 3.
Ahomo = 3;

[x,error] = fp_solver(mesh, At, Ahomo, L);
[cf,cof] = compute_cf(error);

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

% Plot the error 
figure;
title(['Error r =',num2str(r),' Lambda=',num2str(L),' Cont.Factor=',num2str(cf)]);
semilogy(error,'ko--');

% plot contraction factor based error esitmate
b2 = cof(1)*[1:10] + cof(2);
hold on;semilogy(exp(b2),'rd:');

legend('H1-semi. error','error based on contraction factor');
