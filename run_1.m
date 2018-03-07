%
% Antti Hannukainen 4.3.2018 / Otaniemi
%
%----------------------------------------------------------------- 
%
% Solve the Poisson equation (A \nabla u,\nabla v) = (f,v) on domain 
% (0,r)x(0,r) with zero dirichlet boundary condition and load function 
% f = 1. The coefficient function A in each cell has value 1 / 9 
% with 50 / 50 probability. Finite element method with triangular 
% P1-elements is used.
%
% The linear system is solved using fixed point iteration from Theorem 1.1.
% The parameter Lambda is fixed and the size of the domain is varied
% between [10:20:100]. Due to the random nature of the problem an
% averaged contraction factor is computed for each r from ten simulation
% runs. This corresponds to test in Figure 4.
% 
% requires util-folder in the path
% 

close all;
clear all;

r_list = 10:20:100; % values of r used in the computation
Nref = 1;           % number or refinements for the FE - mesh
Nave = 10;          % number of averaging steps used. 
   L = 0.4;         % Lambda parameter for iterative solver.
   
for i=1:length(r_list)
   
    r = r_list(i);
    
    % generate mesh
    [mesh,t2c] = make_Ur_mesh(r,Nref);
    
    for n=1:Nave
    
        % generate random pwc. on each cell of Ur.
        At = make_1_9_cell_At(r, t2c);

        % The corresponding homogenised parameter is 3.
        Ahomo = 3;

        [x,error] = fp_solver(mesh, At, Ahomo, L);
        cf(i,n) = compute_cf(error);

    end
    
    cfave(i) = mean(cf(i,:));
end
  
% Plot the averaged convergence factor
figure;
title(['Averaged contraction factor for L=',num2str(L)]);
plot(r_list,cfave,'ko--');
xlabel('r'); ylabel('averaged contraction factor');
