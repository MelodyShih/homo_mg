% 
% Antti Hannukainen and Mika Juntunen 1.3.2007 
% 
%-----------------------------------------------------------------
%
% affine_tet.m constructs vectorized affine mapping.
%
% The vectorized affine mapping consists out of matrices
% Ax, Ay, Px, Py and vectors bx, by. Row i in Ax(i,:) corresponds
% to affine mapping of component x from the reference element to
% element i. Row i in Px(i,:) corresponds to inverse transpose of 
% the affine mapping matrix of component x from the reference element to 
% element i. 
%
% The affine mapping from reference coordinates (x,y)
% to element i is 
%
%  Ax(i,:) *  x  +  bx(i)
%  Ay(i,:)    y     by(i)
%
% The structure is designed to support operations for all
% elements simultaniously. In practice this means, that a mapping
% of integration points from reference element to global element
% can be done simply as
%
% Ax *  x1 x2 x3 ...   + bx
%       y1 y2 y3 ...
%
% Ay *  x1 x2 x3 ...   + by
%       y1 y2 y3 ..
%
% The action of the inverse of the affine mapping from reference coordinates 
% (x,y) to element i is collected to matrices Px and Py. 
%
% The calling syntax of affine_tri is
%
% function map = affine_tri(mesh)
%
%   mesh    = mesh structure representing the triangulation, see inittri.m
%
%   map.Ax  =  Nt x 2 matrices (see above)
%   map.Ay  =  Nt x 2 matrices (see above)
%
%   mat.bx   = Nt x 1 vector (see above)
%   mat.by   = Nt x 1 vector (see above)
%   
%   map.Px  =  Nt x 2 matrices (see above)
%   map.Py  =  Nt x 2 matrices (see above)
%
%   mat.detA = Nt x 1 vector. detA(i) is the determinant of the matrix A 
%              in affine mapping  to element i. 
%   




function map = affine_tri(mesh)



map.Ax = [(mesh.p(1,mesh.t(2,:))-mesh.p(1,mesh.t(1,:)))' ...
      (mesh.p(1,mesh.t(3,:))-mesh.p(1,mesh.t(1,:)))'];

map.Ay = [(mesh.p(2,mesh.t(2,:))-mesh.p(2,mesh.t(1,:)))' ...
      (mesh.p(2,mesh.t(3,:))-mesh.p(2,mesh.t(1,:)))'];


map.bx = mesh.p(1,mesh.t(1,:))';
map.by = mesh.p(2,mesh.t(1,:))';


% Determinant
%-------------
map.detA = -map.Ax(:,2).*map.Ay(:,1) + map.Ax(:,1).*map.Ay(:,2) ;

% Transpose of inverse of A
%----------------------------
map.Px = bsxfun(@rdivide,[ map.Ay(:,2) -map.Ay(:,1) ],map.detA);
map.Py = bsxfun(@rdivide,[ -map.Ax(:,2)  map.Ax(:,1)],map.detA);







