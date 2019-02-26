%
% Melody Shih 02.18.2018
% (Adapted from run_1 by Antti Hannukainen 4.3.2018 / Otaniemi)
%
%----------------------------------------------------------------- 
% Study the effect of contrast of coefficient, i.e. cmax/cmin
% 
% The homogenized coefficient is fixed to 3.
%
% Solve the Poisson equation (A \nabla u,\nabla v) = (f,v) on domain 
% (0,r)x(0,r) with zero dirichlet boundary condition and load function 
% f = 1. The coefficient function A in each cell has value cmin / cmax
% with 50 / 50 probability. Finite element method with triangular 
% P1-elements is used.
% 
% requires util-folder in the path
% 
% 

close all;
clear;

            r = 32;        % values of r used in the computation
         Nref = 1;          % number or refinements for the FE - mesh
         Nave = 5;          % number of averaging steps used. 
            L = 0.4;        % Lambda parameter for iterative solver.
contrast_list = 10.^(1:4);
         iter = zeros(size(contrast_list));
        chomo = 9;
        Ahomo = [3.164095 4.519 11 36.31054];
            p = 0.6;        % p for cmax

for i=1:length(contrast_list)
      
    contrast = contrast_list(i);
    cmin = sqrt(chomo/contrast);
    cmax = cmin*contrast;
    fprintf('cmax = %f, cmin = %f\n', cmax, cmin);
    
    % generate mesh
    [mesh,t2c] = make_Ur_mesh(r,Nref);
    
    for n=1:Nave
    
        % generate random pwc. on each cell of Ur.
        At = make_cmin_cmax_cell_At(r, t2c, cmin, cmax, p);

        % The corresponding homogenised parameter is 3.
        ahomo = Ahomo(i);

        [x,error,iterone] = fp_solver(mesh, At, ahomo, L);
        iter(i) = iter(i) + iterone;
        cf(i,n) = compute_cf(error);

    end
    iter(i) = iter(i)/Nave;
    cfave(i) = mean(cf(i,:));
end
%% Plots
close all;

% Plot the averaged convergence factor
figure;
semilogx(contrast_list,iter,'ko--');
title(['Number of Iterations to achieve $10^{-9}$ error in H1-semi norm, $L =$',num2str(L), ', $r =$', num2str(r)], ...
      'FontSize',12,'Interpreter','latex');
xlabel('$c_{max}/c_{min}$', 'FontSize',15,'Interpreter','latex'); 
ylabel('Iterations','FontSize',15,'Interpreter','latex');
grid on;

% Plot the averaged convergence factor
figure;
semilogx(contrast_list,cfave,'ko--');
title(['$\lambda=$',num2str(L), ', $r = $', num2str(r)], ...
      'FontSize',12,'Interpreter','latex');
xlabel('$c_{max}/c_{min}$','FontSize',15,'Interpreter','latex'); 
ylabel('$\rho$','FontSize',15,'Interpreter','latex');
grid on;
