function [d, mesh, xx, qq, p, p2, fluxmax] = Project2(h, eltypestr, r)
    % Project 2
    

    %Given values
    TopFlux = 3; %W/mm, vector
    TRight = 0; %C
    thickness = 2; %mm
    Tair = 20; %C
    k = 0.4; %W/mm-K
    convection = 0.0001; %W/mm^2-K

    if nargin == 0
        h = 1;
        eltypestr = 'quad8';
        r = 1; %mm
    end

    path = append('./',eltypestr,'-r',string(r),'-h',string(h),'.inp');

    mesh = abaqus_reader(path);

    d = zeros(length(mesh.x),1);
    f = zeros(length(mesh.x),1);

    D = eye(2)*k;

    K = spalloc(length(mesh.x),length(mesh.x),9*length(mesh.x));

    if strcmp(eltypestr, 'quad4')
        qpts = [[-1 1 1 -1; -1 -1 1 1]/sqrt(3); 1 1 1 1];
    elseif strcmp(eltypestr, 'quad8')
        qpts = zeros(3,9);
        linquad3 = [-sqrt(3/5), 0, sqrt(3/5); 5/9, 8/9, 5/9];
        m = 1;
        for i = 1:3
            for j = 1:3
                ksi = linquad3(1,i);
                eta = linquad3(1,j);
                weight = linquad3(2,i)*linquad3(2,j);
                qpts(:,m)= [ksi;eta;weight];
                m = m + 1;
            end
        end
    elseif strcmp(eltypestr, 'tri3')
        qpts = [[1 4 1; 1 1 4]/6; [1 1 1]/6];
    elseif strcmp(eltypestr, 'tri6')
        qpts = [0.1012865073, 0.7974269853, 0.1012865073, 0.4701420641, 0.4701420641, 0.0597158717, 0.3333333333;
                0.1012865073, 0.1012865073, 0.7974269853, 0.0597158717, 0.4701420641, 0.4701420641, 0.3333333333;
                0.0629695903, 0.0629695903, 0.0629695903, 0.0661970764, 0.0661970764, 0.0661970764, 0.1125000000];
    end
        
    for c = mesh.conn

        xe = mesh.x(:,c);
        Ke = zeros(length(c));
        fe = zeros(length(c),1);
        Kh = zeros(length(c));

        for q = qpts

            if strcmp(eltypestr, 'quad4')
                [N, dNdp] = shape_quad4(q);
            elseif strcmp(eltypestr, 'quad8')
                [N, dNdp] = shape_quad8(q);
            elseif strcmp(eltypestr, 'tri3')
                [N, dNdp] = shape_tri3(q);
            elseif strcmp(eltypestr, 'tri6')
                [N, dNdp] = shape_tri6(q);
            end

            J = xe*dNdp;
            dNdx = dNdp/J;
            Ke = Ke + dNdx*D*dNdx'*det(J)*q(end)*thickness;
            fe = fe + N*convection*Tair*det(J)*q(end);
            Kh = Kh + N*convection*N'*det(J)*q(end);

        end

        K(c,c) = K(c,c) + Ke + Kh;
        f(c) = f(c) + fe;

    end

    ns_top = find(mesh.x(2,:) == max(mesh.x(2,:)));
    [~,order] = sort(mesh.x(1,ns_top));
    ns_top = ns_top(order);
    
    if strcmp(eltypestr,'tri3') || strcmp(eltypestr,'quad4')
       fluxset = [ns_top(1:end-1);
                  ns_top(2:end)];
    else
       fluxset = [ns_top(1:2:end-2);
                  ns_top(2:2:end-1);
                  ns_top(3:2:end)];
    end
    
    for c = fluxset

        xe = mesh.x(:,c);
        
        quadpts = [[-1 1]/sqrt(3);[1 1]];
        for q = quadpts
            
            if strcmp(eltypestr,'tri3') || strcmp(eltypestr,'quad4')
                [N, dNdp] = shape2(q(1));
            else 
                [N, dNdp] = shape3(q(1));
            end
            J = xe*dNdp;
            f(c) = f(c) + N*TopFlux*norm(J)*q(end);
        end


    end

    ns_right = find(mesh.x(1,:) == max(mesh.x(1,:)));
    K(ns_right,:) = 0;
    K(ns_right, ns_right) = eye(length(ns_right));
    f(ns_right) = TRight;

    d = K\f;
    
    
    A = spalloc(length(mesh.x),length(mesh.x),9*length(mesh.x));
    y = zeros(length(mesh.x),2);
    for c = mesh.conn
        xe = mesh.x(:,c);
        de = d(c)';
        Ae = zeros(length(c));
        for q = qpts
            if strcmp(eltypestr,'quad4')
                [N, dNdp] = shape_quad4(q);
            elseif strcmp(eltypestr,'quad8')
                [N, dNdp] = shape_quad8(q); 
            elseif strcmp(eltypestr,'tri3')
                [N, dNdp] = shape_tri3(q);
            elseif strcmp(eltypestr,'tri6')
                [N, dNdp] = shape_tri6(q); 
            end
            J = xe*dNdp;
            dNdx = dNdp/J;
            qh = -k*de*dNdx;
            Ae = Ae + N*N'*det(J)*q(end);
            y(c,:) = y(c,:) + N*qh*det(J)*q(end);
        end
        A(c,c) = A(c,c) + Ae;
    end
    qint = A\y;
    q = (qint(:,1).^2+qint(:,2).^2).^0.5;

    
    
    p.vertices = mesh.x';
    p.faces = mesh.conn';
    p.facecolor = 'interp';
    p.facevertexcdata = d;
    
    if strcmp(eltypestr, 'quad8')
        p.faces = mesh.conn([1;5;2;6;3;7;4;8],:)';
    elseif strcmp(eltypestr, 'tri6')
        p.faces = mesh.conn([1;4;2;5;3;6],:)';
    end
    
    p2.vertices = mesh.x';
    p2.faces = mesh.conn';
    p2.facecolor = 'interp';
    p2.facevertexcdata = q;
    fluxmax = max(q);
    
    if strcmp(eltypestr, 'quad8')
        p2.faces = mesh.conn([1;5;2;6;3;7;4;8],:)';
    elseif strcmp(eltypestr, 'tri6')
        p2.faces = mesh.conn([1;4;2;5;3;6],:)';
    end
    
    xx = [];
    qq = [];
    
    for c = mesh.conn
       
        xe = mesh.x(:,c);
        de = d(c)';
        if strcmp(eltypestr,'quad4')
           [N, dNdp] = shape_quad4([0;0]);
        elseif strcmp(eltypestr,'quad8')
           [N, dNdp] = shape_quad8([0;0]); 
        elseif strcmp(eltypestr,'tri3')
           [N, dNdp] = shape_tri3([1/3;1/3]);
        elseif strcmp(eltypestr,'tri6')
           [N, dNdp] = shape_tri6([1/3;1/3]); 
        end
        J = xe*dNdp;
        dNdx = dNdp/J;
        qq(end+1,:) = -(de*dNdx)*D;
        xx(end+1,:) = xe*N;
        
    end

end

% Reads the ABAQUS input file.
% Returns a mesh struct with fields
% mesh.x   : NSDxNN array of nodal coordinates
% mesh.nid : global node ids (usually 1:NN)
% mesh.nn  : number of nodes in the mesh
% mesh.conn: connectivity matrix NNExNE
% mesh.eid : global element ids (not always 1:NE)
function [mesh] = abaqus_reader(path)
    fid = fopen(path);
    mesh.x = [];
    mesh.nid = [];
    mesh.conn = [];
    mesh.eid = [];

    % Information required for sidesets (boundary information).
    mesh.elset = {};
    mesh.elset_id = {};
    mesh.surface = {};
    mesh.surface_id = {};

    line = fgetl(fid);
	while ~feof(fid)
        if strfind(line, '*NODE, NSET=ALLNODES') == 1
            % Reads nodal data.
            while ~feof(fid)
                line = fgetl(fid); 
                if strcmp(line,'**')
                    break
                else
                    % Can read either 2 or 3 floats.
                    [x, ct] = sscanf(line, '%d, %f, %f, %f');
                    assert(ct==3 || ct==4, 'Error reading mesh');
                    mesh.nid(end+1) = x(1);
                    mesh.x(:,end+1) = x(2:end);
                end
            end
            line = fgetl(fid);
        % Reads connectivity data.
        elseif strfind(line, '*ELEMENT, TYPE=')
            while ~feof(fid)
                line = fgetl(fid); 
                if strcmp(line,'**')
                    break
                else
                    [x, ct] = sscanf(line, '%d,');
                    assert(ct>3, 'Error reading mesh');                    
                    mesh.eid(end+1)    = x(1);
                    mesh.conn(:,end+1) = x(2:end);
                end
            end
            line = fgetl(fid);
        % Reads element sets and stores their elements and names in cell arrays.
        elseif strfind(line, '*ELSET, ELSET=')
            line = strsplit(line, '=');            
            mesh.elset_id{end+1} = line{end};
            mesh.elset{end+1} = [];
            while ~feof(fid)
                line = fgetl(fid);
                if ~isempty(line) && line(1) == '*'
                    break
                else                    
                    [x,ct] = sscanf(line, '%d,');
                    mesh.elset{end} = [mesh.elset{end}; x];
                end
            end
        % Reads surfaces information (sidesets from cubit).
        elseif strfind(line, '*SURFACE, NAME=')
            line = strsplit(line, '=');            
            mesh.surface{end+1} = {};
            mesh.surface_id{end+1} = line{end};
            while ~feof(fid)
                line = fgetl(fid);
                if ~isempty(line) && line(1) == '*'
                    break
                elseif ~isempty(line)
                    % Line should look like [elset_id], [edge number]
                    line = strsplit(line, ',');
                    assert(length(line)==2, 'Bad surface information')
                    [x,ct] = sscanf(line{2}, ' E%d');
                    mesh.surface{end}{end+1} = {line{1}, x};
                end
            end
        else
            line = fgetl(fid);
        end        
    end
    mesh.nn = size(mesh.x,2);
    mesh.ne = size(mesh.conn,2);
    
    % Correction for 4-point tet meshes.
    if size(mesh.conn, 1) == 4 && size(mesh.x, 1) == 3
        mesh.conn = mesh.conn([2,1,3,4],:);
    end

    mesh.edge_connectivity = @(id) build_edge_connectivity(mesh, id);
    mesh.get_surface_faces = @(id) build_faces(mesh, id);
    mesh.closest_node      = @(pt) closest_node(mesh, pt);

    fclose(fid);
    fprintf('Read %d nodes and %d elements from %s.\n',...
        mesh.nn, mesh.ne, path);
end

% Returns an m by n array giving the connectivity matrix of edges along
% named element sidesets, where m is the number of node along an element 
% edge and n is the number of element edges on the sideset.  This is used
% for defining heat flux or traction boundary conditions on the boundary
% of the domain.
function [conn] = build_edge_connectivity(mesh, id)
    if size(mesh.conn, 1) == 3
        edge_nodes = [1 2 3; 2 3 1];
    elseif size(mesh.conn, 1) == 4
        edge_nodes = [1 2 3 4; 2 3 4 1];
    elseif size(mesh.conn, 1) == 6
        edge_nodes = [1 2 3; 4 5 6; 2 3 1];
    else
        error('abaqus_reader: cannot read mesh type.');
    end
    conn = [];
    for i = 1:length(mesh.elset_id)        
        pattern = strcat(id, '_S(\d)');
        match = regexp(mesh.elset_id{i}, pattern, 'tokens');
                
        if ~isempty(match)
            edge = str2num(match{:}{1});
            conn = [conn, mesh.conn(edge_nodes(:, edge), mesh.elset{i})];
        end
    end
end

% Returns the index of the nearest node to pt. (somewhat slow).
function [n] = closest_node(mesh, pt)
    dx = mesh.x(1,:) - pt(1);
    dy = mesh.x(2,:) - pt(2);    
    [rmin, n] = min(dx.*dx + dy.*dy);
end

function [N, dNdp] = shape_quad4(p)
%Returns the shape functions and their gradients for a 4-node element
N = (1/4)*[(1-p(1))*(1-p(2));
           (1+p(1))*(1-p(2));
           (1+p(1))*(1+p(2));
           (1-p(1))*(1+p(2))]; % (4x1) replace with correct expression
dNdp = (1/4)*[ -1+p(2), -1+p(1);
               1-p(2), -1-p(1);
               1+p(2), 1+p(1);
               -1-p(2), 1-p(1)]; % (4x2) replace with correct expression
end

function [N, dNdp] = shape_quad8(p)
%Returns the shape functions and their gradients for an 8-node element
N = (1/4)*[-(1-p(1))*(1-p(2))*(1+p(1)+p(2));
           -(1+p(1))*(1-p(2))*(1-p(1)+p(2));
           -(1+p(1))*(1+p(2))*(1-p(1)-p(2));
           -(1-p(1))*(1+p(2))*(1+p(1)-p(2));
           2*(1-p(1)^2)*(1-p(2)); 
           2*(p(1)+1)*(1-p(2)^2);
           2*(1-p(1)^2)*(1+p(2));
           2*(1-p(1))*(1-p(2)^2)]; % (8x1) replace with correct expression
dNdp = (1/4)*[ 2*p(1)+p(2)-2*p(1)*p(2)-p(2)^2, 2*p(2)+p(1)-p(1)^2-2*p(1)*p(2);
               2*p(1)-p(2)-2*p(1)*p(2)+p(2)^2, 2*p(2)-p(1)-p(1)^2+2*p(1)*p(2);
               2*p(1)+p(2)+2*p(1)*p(2)+p(2)^2, 2*p(2)+p(1)+p(1)^2+2*p(1)*p(2);
               2*p(1)-p(2)+2*p(1)*p(2)-p(2)^2, 2*p(2)-p(1)+p(1)^2-2*p(1)*p(2);
               4*p(1)*p(2)-4*p(1), 2*p(1)^2-2;
               -2*p(2)^2+2, -4*p(1)*p(2)-4*p(2);
               -4*p(1)*p(2)-4*p(1), -2*p(1)^2+2;
               2*p(2)^2-2, 4*p(1)*p(2)-4*p(2)]; % (8x2) replace with correct expression
end

function [N, dNdp] = shape_tri3(p)
%Returns the shape function and their gradients for a 3-node element.
N = [ p(1);
      p(2);
      1-p(1)-p(2)]; %Replace with correct expression
dNdp = [ 1, 0;
         0, 1;
        -1, -1]; %Replace with correct expression
end

function [N, dNdp] = shape_tri6(p)
%Returns the shape function and their gradients for a 6-node element.
N = [ p(1)*(2*p(1)-1);
      p(2)*(2*p(2)-1);
      (1-p(1)-p(2))*(2*(1-p(1)-p(2))-1);
      4*p(1)*p(2);
      4*p(2)*(1-p(1)-p(2));
      4*p(1)*(1-p(1)-p(2))]; %Replace with correct expression
dNdp = [ 4*p(1) - 1, 0;
         0, 4*p(2) - 1;
         4*p(1)-3+4*p(2), -3+4*p(1)+4*p(2);
         4*p(2), 4*p(1);
         -4*p(2), 4-4*p(1)-8*p(2);
         4-8*p(1)-4*p(2), -4*p(1)]; %Replace with correct expression
end

function [N, dNdp] = shape3(p)
 N = [(-1/2)*p*(1-p) ; (1-p)*(1+p) ; (1/2)*p*(1+p)];
 dNdp = [p - 0.5 ; -2*p ; p + 0.5];
end

function [N, dNdp] = shape2(p)
 N = (1/2)*[1-p;1+p];
 dNdp = [-0.5;0.5];
end