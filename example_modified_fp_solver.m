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
% clear all;

      r = 10;  % size of the domain Ur = (0,r)x(0,r)
   Nref = 1;   % number or refinements for the FE - mesh
   cmin = 1;
   cmax = 9;
      w = 2/3; % weight for weighted jacobian smoother
maxiter = 50;  
   
% generate mesh
[mesh,t2c] = make_Ur_mesh(r,Nref);

% generate random pwc. on each cell of Ur.
At = make_cmin_cmax_cell_At(r, t2c, cmin, cmax);

% The corresponding homogenised parameter is 3.
Ahomo = sqrt(cmin*cmax);



[x,error,iter,~,~] = modified_fp_solver(mesh, At, Ahomo, 0, w, maxiter, 1);
% [x1,error1,iter1,~,~] = jacobi_solver(mesh, At, Ahomo, 50, w);

[cf,cof] = compute_cf(error);
% [cf1,cof1] = compute_cf(error1);

%% Plot Result
close all;

% Plot the solution.
X = mesh.p(1,:);
Y = mesh.p(2,:);
t = mesh.t;

% figure;
title(['Solution $u, r =$',num2str(r)],'Fontsize',15,'Interpreter','latex');
P = patch(X(t),Y(t),x(t),x(t));
set(P,'edgecolor','none');
hold on;
colorbar;
axis equal;

% Plot the coefficient field
figure;
title(['Coefficient field, r =',num2str(r)],'Fontsize',15,'Interpreter','latex');
P = patch(X(t),Y(t),repmat(At(:)',3,1),repmat(At(:)',3,1));
set(P,'edgecolor','none');
xlabel('$x$', 'Fontsize',15,'Interpreter','latex');
ylabel('$y$', 'Fontsize',15,'Interpreter','latex');
axis square;

% Plot the error 
figure;
semilogy(error,'ko-');
hold on;
semilogy(error1,'bo-')
title(['H1-semi. Error, $r =$',num2str(r)], ...
    'Fontsize',15,'Interpreter','latex');
xlabel('Iterations', 'Fontsize',15,'Interpreter','latex');
lgd = legend('Homogenized Solve + Jacobian smoother', 'Jacobi Iterative');
lgd.Interpreter = 'latex';
lgd.FontSize = 10;
lgd.Location = 'Best';
grid on;