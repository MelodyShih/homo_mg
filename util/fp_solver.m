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


function [xFP,error] = fp_solver(mesh, At, Ahom, Lmassive)

disp('assembly');

% The exact problem.
[A,b] = assembly_P1(mesh,At,0,@(x,y)(ones(size(x))));

% The homogenised problem
AH = assembly_P1(mesh,Ahom,0,@(x,y)(ones(size(x))));

% The exact problem with large L2-term.
AM =  assembly_P1(mesh,At,Lmassive^2,@(x,y)(ones(size(x))));

% The homogenised problem with large L2-term.
AHM =  assembly_P1(mesh,Ahom,Lmassive^2,@(x,y)(ones(size(x))));

 
 M = assembly_P1(mesh,0,1,@(x,y)(ones(size(x))));
A0 = assembly_P1(mesh,1,0,@(x,y)(ones(size(x))));

disp('assembly done');

% Zero initial guess
xFP = 0*b;
 in = mesh.in;

disp('exact solve');
xE = 0*b;
xE(in) = A(in,in)\b(in);
disp('exact solve done');

% evaluate the for initial guess xFP = 0.
error(1) = sqrt( (xE)'*A0*(xE) );
 
% precompute Cholesky factorisations.    
disp('factorisations (rev)');

[RM,p,SM] = chol(AM(in,in));
RMT = RM';

[RH,p,SH] = chol(AH(in,in));
RHT = RH';

disp('factorisations done (rev)');

for i=1:10
     
    %% Compute u_0 : 
    res1 = b - A*xFP;
    u0 = 0*b;
    u0(in) = SM*(RM\(RMT\(SM'*res1(in))));
       
    %% Compute u_hom
    res2 = Lmassive^2*M*u0;
    
    uhom = 0*b;
    uhom(in) = SH*(RH\(RHT\(SH'*res2(in))));

    %% utilde
    res3 = AHM*uhom;
    utilde =  0*b;
    utilde(in) = SM*(RM\(RMT\(SM'*res3(in))));

  
    %% Update
    xFP(in) = xFP(in) + u0(in)  + utilde(in);
        
    % evaluate error 
    error(i+1) = sqrt( (xFP-xE)'*A0*(xFP-xE) );
    
    disp(['fixed-point iteration step:',num2str(i),'error (H1-semi) :',num2str(error(i+1))]);
        
end
