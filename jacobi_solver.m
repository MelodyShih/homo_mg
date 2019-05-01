%
% Antti Hannukainen 4.3.2018 / Otaniemi
%
% Revision 17.4.2018 / Otaniemi (Simplify the method)
%
%----------------------------------------------------------------- 
%
% Implementation of the fixed point iteration described in Th 1.1
% 
% This version is meant to study the behaviour of the contraction factor and
% uses exact inverses. Ten steps of the iteration are used starting from zero initial guess
% The load is set to one.
%
% CALLING SYNTAX IS
%
% function [xFP,error] = fp_solver(mesh, At, Ahomo, Lmassive)
% 
%   mesh    = mesh structure representing the triangulation, see inittri.m
%  
%   At      = Nt-vector of the value of coefficient A over mesh. 
%             A_|K = At(K). Scalar At is appropriately extended to vector of
%             correct size.
%
%   Ahomo   = Homogenised coefficient correspondig to At. 
% 
%   Lmassive = The massive term.
%
%       xFP = The solution of the fixed-point iteration after ten steps.
%
%     error = 10x1 - vector of the H1-seminorm of the error between exact 
%             FE-solution and each of the iterates.
%


function [xFP,error,i,R,Din] = jacobi_solver(mesh, At, Ahom, nsweeps, w)

disp('assembly');

% The exact problem.
[A,b] = assembly_P1(mesh,At,0,@(x,y)(ones(size(x))));


% The homogenised problem
% AH = assembly_P1(mesh,Ahom,0,@(x,y)(ones(size(x))));

% The exact problem with large L2-term.
% AM =  assembly_P1(mesh,At,Lmassive^2,@(x,y)(ones(size(x))));

% The homogenised problem with large L2-term.
% AHM =  assembly_P1(mesh,Ahom,Lmassive^2,@(x,y)(ones(size(x))));

 
%  M = assembly_P1(mesh,0,1,@(x,y)(ones(size(x))));
A0 = assembly_P1(mesh,1,0,@(x,y)(ones(size(x))));

disp('assembly done');

% Zero initial guess
in = mesh.in;

% Set the right hand side such that the true solution has a_i = 1
% [E,~] = eig(full(A(in,in)));
% [~,n] = size(E);
% reductionfac = zeros(1,n);
% b = 0*b;

% Initial guess
xFP = 0*b;
% initialize x such that the error in each eigendirection is 1.
% xFP(in) = E*ones(size(b(in)));

Din = diag(diag(A(in,in)));
R = A(in,in) - Din;

disp('exact solve');
xE = 0*b;
xE(in) = A(in,in)\b(in);
disp('exact solve done');

maxiter = 50;
error = zeros(1,maxiter);

% evaluate the for initial guess xFP = 0.
error(1) = sqrt( (xE)'*A0*(xE) );
 
% precompute Cholesky factorisations.    
% disp('factorisations (rev)');

% [RM,p,SM] = chol(AM(in,in));
% RMT = RM';

% [RH,~,SH] = chol(AH(in,in));
% RHT = RH';

% disp('factorisations done (rev)');

% totaliter = 0;
disp(['fixed-point iteration step:',num2str(0),'error (H1-semi) :',num2str(error(1))]);
for i=1:maxiter   
    %% Pre-smoothing: 
    for s=1:nsweeps
        xFP(in) = w*(Din\(b(in) - R*xFP(in))) + (1-w)*xFP(in);
    end
%     eigerrcomp = E\xFP(in);
%     reductionfac = abs(eigerrcomp).^(1/nsweeps);
%     figure;
%     plot(reductionfac, '-o')
       
    %% Compute u_hom
%     res2 = b(in) - A(in,in)*xFP(in);
%     uhom = 0*b;
%     uhom(in) = SH*(RH\(RHT\(SH'*res2)));
    
    %% Update x
%     xFP(in) = xFP(in)+uhom(in);
    
    %% Post-smoothing
%     for s=1:nsweeps
%        xFP(in) = w*(Din\(b(in) - R*xFP(in))) + (1-w)*xFP(in);
%     end
    
    %% evaluate error 
    error(i+1) = sqrt( (xFP-xE)'*A0*(xFP-xE) );
    disp(['fixed-point iteration step:',num2str(i),'error (H1-semi) :',num2str(error(i+1))]);
    if(error(i+1) < 1e-10)
        break;
    end    
end
end

