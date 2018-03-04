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
% The size of the domain is fixed and parameter Lambda is varied
% logarihmically between values 0.01 and 0.5. Due to the random nature 
% of the problem an averaged conraction factor is computed for each 
% Lambda. This corresponds to test in Figure 4.
% 

clear all;
close all;

r = 100;      % value  of r used in the computation
Nref = 1;     % number or refinements for the FE - mesh
Nave = 2;     % number of averaging steps used. 
L_list = logspace(log10(0.01),log10(0.5),10); % Lambda parameter.

% generate mesh
[mesh,t2c] = make_Ur_mesh(r,Nref);  

for i=1:length(L_list)
    
    L = L_list(i);
    
    for n=1:Nave
    
        % generate random pwc. on each cell of Ur.
        At = make_1_9_cell_At(r, t2c);

        % The corresponding homogenised parameter is 3.
        Ahomo = 3;

        [x,error] = fp_solver(mesh, At, Ahomo, L);
        cf(i,n) = compute_cf(error);

    end
    
    % averaged contraction factor.
    cfave(i) = mean(cf(i,:));
end
  
% Plot the averaged convergence factor
figure;
title(['exponential of the av. contraction factor for r=',num2str(r)]);
xlabel('\lambda');
ylabel('exp(\rho)');
loglog(L_list,exp(cfave),'ko--');
hold on;
loglog(L_list,sqrt(L_list),'k--');
legend('computed','\lambda^{1/2}','Location','southeast');
