% Antti Hannukainen / 11.1.2008 / Otaniemi
%
% plot_2Dtri_mesh.m
%
% Plot 2D triangular mesh. The plotting is done by using line or patch -
% functions. Using patch function is faster, but might result into low 
% quality. For more information on mesh - data structure, see inittri.m 
%
% The calling syntax of plot_2Dtri_mesh.m is
%
% H = plot_2Dtri_mesh(mesh,method) 
%
% mesh         =  mesh-structure to be plotted (must be 2D triangles)
%
% method       =  method can be choosen between 'lines' / 'patch' (default).
%                 Using option 'patch' can produce low quality .eps - files, 
%                 but is much faster. The 'line' option is recommended when 
%                 the aim is to export .eps - files
%
% RETURNS
% 
% H            = Handle to drawn graphic object ( can be used to change
%                                                 color etc. see set(H) )
%

function H=plot_2Dtri_mesh(mesh,method) 

if(nargin ==1)
    method = '';
end

if( strcmp(method,'lines'))

    X(1,:) = mesh.p(1,mesh.edges(1,:));
    X(2,:) = mesh.p(1,mesh.edges(2,:));
    
    Y(1,:) = mesh.p(2,mesh.edges(1,:));
    Y(2,:) = mesh.p(2,mesh.edges(2,:));
   
    H = line(X,Y);
    set(H,'Color',[0 0 0]);
else
    
    px = mesh.p(1,:);
    py = mesh.p(2,:);
    
    X = px(mesh.t);
    Y = py(mesh.t);
    
    H= patch(X,Y,0*X);set(H,'FaceColor','none');

end


