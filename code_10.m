function [f] = code_10(mesh, qbar)
    if nargin == 0
        mesh = abaqus_reader('semicircle.inp');
        qbar = @(x,y) x.^2 - y.^2;
    end
end

% Finds nodes located at radius of 10 or larger from the origin.
ns = find(sqrt(mesh.x(1,:).^2 + mesh.x(2,:).^2) > 9.999);
% Reorders the nodes so that they are sorted by their x-coordinate.
[~, order] = sort(mesh.x(1,ns));
ns = ns(order);


% Code used to generate figure (Not needed for submission)
r = sqrt(mesh.x(1,:).^2 + mesh.x(2,:).^2);
ns = find(r > 9.999);
figure(1); clf; hold on;
p.facecolor = 'w';
p.vertices = mesh.x';
p.faces = mesh.conn([1 5 2 6 3 7 4 8], :)';
patch(p);
scatter(mesh.x(1,ns), mesh.x(2,ns));
    for i=ns
        text(mesh.x(1,i), mesh.x(2,i), int2str(i));
    end 
axis equal;
