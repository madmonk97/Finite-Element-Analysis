
function [d, sig] = code_3(L, H, n, p)

    % Defining the Input arguments - L H n p
    if nargin == 0
        clc
        L = 5;      % Length of Truss (m)
        H = 1;      % Height of Truss (m)
        n = 5;      % Nuber of segments along truss (n=5)(n>2)
        p = 10;     % magnitude of the three concentrated forces (N)
    end
    
    mesh = make_truss(L, H, n);
    
    % Specify the given data values 
    E   = 210e9;                        % in Pa
    A   = pi/4*((0.02)^2-(0.012)^2);    % in m^2
    g   = 9.81;                         % Gravity in m/s^2
    rho = 8050;                         % Density in kg/m^3           
    
    K  = zeros(2*length(mesh.x));      % Create space for Stiffness matrix
    f  = zeros(2*length(mesh.x),1);    % Create space for Force matrix
    
    % Addition of external forces to the given 3 rightmost bottom points in the truss figure
    f(length(f)/2)  = f(length(f)/2) - p;
    f(length(f)/2-2)= f(length(f)/2-2) - p;
    f(length(f)/2-4)= f(length(f)/2-4) - p;
       
    for c           = mesh.conn
        xe          = mesh.x(:,c);                       % Nodal position
        dx          = xe(:,2) - xe(:,1);                 % Difference between 2 nodes (node 2 - node 1)
        t           = dx/norm(dx);                       % Dividing dx by norm(dx) (absolute values)
        Re          = [t' 0 0; 0 0 t'];                  % Rotation matrix
        Ke          = Re'*E*A/norm(dx)*[1 -1;-1 1]*Re;   % Element stiffness matrix
        sctr(1:2:4) = 2*c-1;                             % x- degree of freedom is 2*c-1
        sctr(2:2:4) = 2*c;                               % y- degree of freedom is 2*c
        w           = rho*g*norm(dx)*A;                  % weight = density*g*length of element*Area
        f(2*c)      = f(2*c)- (1/2)*w;                   % Forces on y f(2*c) = f(2*c) + averaged weight
        K(sctr,sctr)= K(sctr,sctr) + Ke;                 % K matrix (Global)
    end
    
    % Defining the boundary conditions
    fixed_dof               = [1*2-1,1*2,(n+2)*2-1]; % Fixing the node 1 in x,y direction and the node 7 in x direction.
    K(fixed_dof,:)          = 0;
    K(fixed_dof,fixed_dof)  = eye(length(fixed_dof));
    f(fixed_dof)            = 0;
    
    d = K\f; % Displacement calculations
    
    sig     = [];
    counter = 1;
    
    for c                   = mesh.conn
        xe                  = mesh.x(:,c);      % Nodal position
        dx                  = xe(:,2) - xe(:,1);% Difference between 2 nodes (node 2 - node 1)
        t                   = dx/norm(dx);      % Dividing dx by norm(dx) (absolute values)
        Re                  = [t' 0 0; 0 0 t']; % Rotation matrix
        sctr(1:2:4)         = 2*c-1;            % x- degree of freedom is 2*c-1
        sctr(2:2:4)         = 2*c;              % y- degree of freedom is 2*c
        
        original_length     = norm(dx);                                     % Original length
        de                  = d(sctr);                                      % Gather - element displacement
        intermediate_step   = Re*de;         
        delta_L             = intermediate_step(2)-intermediate_step(1);    % Change in length (delta)
        epsilon             = delta_L/original_length;                      % Strain
        
        sig(counter) = E*epsilon;
        
        counter = counter + 1;
        
    end
    
end

