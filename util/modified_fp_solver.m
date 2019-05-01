%
% Melody Shih 26.2.2019
% (Adapted from fp_solver by Antti Hannukainen 4.3.2018 / Otaniemi)
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


function [xFP,error,i,R,Din] = modified_fp_solver(mesh, At, Ahom, nsweeps, ...
                                                  w, maxiter, spectral)
if nargin == 4
    maxiter = 30;
    spectral = 0;
end

disp('assembly');

% The exact problem.
[A,b] = assembly_P1(mesh,At,0,@(x,y)(ones(size(x))));


% The homogenised problem
AH = assembly_P1(mesh,Ahom,0,@(x,y)(ones(size(x))));

A0 = assembly_P1(mesh,1,0,@(x,y)(ones(size(x))));

disp('assembly done');

in = mesh.in;

if(spectral)
    [E,V] = eig(full(A(in,in)));
    b = 0*b;
    figure;
    plot(diag(V));
end

% Initial guess
xFP = 0*b;
if(spectral)
    % initialize x such that the error in each eigendirection is 1.
    temp = ones(size(b(in)));
%     temp(10:end) = 0;
    xFP(in) = E*temp;
%     xFP(in) = E*ones(size(b(in)));

end

Din = diag(diag(A(in,in)));
R = A(in,in) - Din;

disp('exact solve');
xE = 0*b;
xE(in) = A(in,in)\b(in);

disp('exact solve done');

error = zeros(1,maxiter);

% evaluate the for initial guess xFP = 0.
error(1) = sqrt( (xE)'*A0*(xE) );
 
% precompute Cholesky factorisations.    
disp('factorisations (rev)');

[RH,~,SH] = chol(AH(in,in));
RHT = RH';

disp('factorisations done (rev)');

% totaliter = 0;
disp(['fixed-point iteration step:',num2str(0),'error (H1-semi) :',num2str(error(1))]);

err1=sqrt( (xFP-xE)'*A0*(xFP-xE) );
for i=1:maxiter   
    %% Pre-smoothing: 
    for s=1:nsweeps
        xFP(in) = w*(Din\(b(in) - R*xFP(in))) + (1-w)*xFP(in);
        if(spectral && s==nsweeps)
            err1 = sqrt( (xFP-xE)'*A0*(xFP-xE) );
%             semilogy([0 length(A(in,in))], [1 1], 'Color','black', 'LineWidth', 2);
%             hold on;
%             title(['Sweeps ', num2str(s),', Err = ',num2str(err1)],'Interpreter', 'latex');
%             temp = E\(xFP(in) - xE(in));   
%             semilogy(abs(temp),'.');
%             ylim([1e-10 100]);
%             grid on;
%             www = waitforbuttonpress;
%             hold off;
        end
    end
    if(spectral)
%         figure;
%         semilogy([0 length(A(in,in))], [1 1], 'Color','black', 'LineWidth', 2);
%         hold on;
%         temp = E\(xFP(in) - xE(in));    
%         semilogy(abs(temp),'.');
%         ylim([1e-16 100]);
    end
    
    %% Compute u_hom
    res2 = b(in) - A(in,in)*xFP(in);
    uhom = 0*b;
    uhom(in) = SH*(RH\(RHT\(SH'*res2)));
    if(spectral)
%         temp = E\uhom(in);
%         semilogy(abs(temp),'x');
    end

    %% Update x
    xFP(in) = xFP(in)+uhom(in);
    if(spectral)
        temp = E\(xFP(in) - xE(in));
        semilogy(abs(temp),'.');
    end
    
    if(spectral)
        lgd=legend('y = 1', 'Pre-smoothing', 'Homogenized Update');%, 'Post-smoothing');
%         lgd=legend('y = 1', 'Homogenized Update');%, 'Post-smoothing');
        lgd.Interpreter='latex';
        lgd.FontSize=10;
        lgd.Location='best';
        err = sqrt( (xFP-xE)'*A0*(xFP-xE) );
        title(['Iteration ', num2str(i),', w = ',num2str(w),', Err=',num2str(err),' (',num2str(err1),')'] ... 
               ,'Interpreter', 'latex');
        grid on;
%         saveas(gcf,['jsmooth', num2str(i),'w',num2str(w*10)],'epsc')
%         saveas(gcf,['homosol'],'epsc');
        ww = waitforbuttonpress;
        hold off;
    end
    
    %% evaluate error 
    error(i+1) = sqrt( (xFP-xE)'*A0*(xFP-xE) );
    disp(['fixed-point iteration step:',num2str(i),'error (H1-semi) :',num2str(error(i+1))]);
    if(error(i+1) < 1e-10)
        break;
    end    
end
end

