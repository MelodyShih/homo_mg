%----------------------------------------------------------------- 
%
% Antti Hannukainen and Mika Juntunen 3.11.2005
%
% Vectorization by Antti Hannukainen 29.8.2007
%      Revision by Antti Hannukainen 27.10.2008
%
%-----------------------------------------------------------------
%
% Construct mesh - structure from a given triangulation. 
%
% The triangulation is expressed using matrices p and t. 
% 
% The p-matrix has dimension 2 x Np. Each column of p represents 
% coordinates of a single node of the mesh. For example p(:,i) 
% is coordinate of node i. Coordinates are assumed to be unique. 
% 
% The t-matrix has dimension 3xNt. Each column of t represents nodal 
% indeces of single triangle in the mesh. For example, t(:,i) is a vector [n1 n2 n3]' 
% such that p(:,n1), p(:,n2), and p(:,n3) are coordinates of the nodes of triangle i. 
%
% MESH STRUCTURE :
% 
%   mesh.p       = nodes in a 2xN-matrix
%
%   mesh.t       = triangles in a 3xN- matrix. 
% 
%   mesh.edges   = a matrix of all edges in the mesh. Each column is an edge 
%                  [n1 ; n2] where n1 and n2 are nodal indeces to p-matrix such that 
%                  n1 < n2; 
%
%   mesh.t2e     = a matrix connecting  triangle's and edges's. 
%                  Each column corresponds to single triangle in the mesh and 
%                  has triangle's edge indeces in the order n1->n2, n2->n3,
%                  n1->n3. 
%
%   mesh.e2t     = inverse of t2e. Boundary edges satisfy mesh.e2t(2,:) = 0.
%
%   mesh.bn      = list of all boundary nodes.
%   mesh.in      = list of all interor nodes.
%
%
% CALLING SYNTAX IS 
%
% function mesh = inittri(p,t)
%
%   p       = points matrix (described above)
%   t       = triangles matrix (described above)
%     
%   mesh    = trimesh structure corresponding to (p,t)
%


function mesh = inittri(p,t)

mesh.p  = p;
mesh.t  = t;

% Initalize size variables
Nt = size(t,2);

e = [1 2; 2 3; 1 3]';

edges = [];

for i=1:size(e,2)
  edges = [edges [ sort( [t(e(1,i),:); t(e(2,i),:)],1)]];
end

[mesh.edges,~,mesh.t2e] = sc_unique(edges);
mesh.t2e = reshape(mesh.t2e,Nt,3)';

% mesh.e2t
e = [mesh.t2e(1,:) mesh.t2e(2,:)  mesh.t2e(3,:)];
t = repmat([1:Nt],1,3);

[ef,If]= unique(e,'first');
[el,Il]= unique(e,'last');

mesh.e2t(1,ef) = t(If);
mesh.e2t(2,el) = t(Il);

mesh.e2t(2,find( (mesh.e2t(1,:)-mesh.e2t(2,:))==0))=0;

be = find(mesh.e2t(2,:) == 0);
mesh.bn = unique([mesh.edges(1,be) mesh.edges(2,be)]);
mesh.in = setdiff(1:size(mesh.p,2),mesh.bn);



