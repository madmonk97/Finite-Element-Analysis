function [d, sigma] = code_4(N, L, w, p)

    % Defining the Input arguments - N, L, w, p
    
    if nargin == 0
        clc
        L = 12;     % length of the truss along the x-direction (m)
        N = 4;      % Height of Truss (m)
        w = 2;      % width of the bottom of the truss
        p = 500;    % Forces/Loads 'p'
    end
    
    mesh = make_3dtruss(N, L, w);
    
    % Specifying the given data in the problem 
    
    E   = 3e9;                          % in Pa
    A   = pi/4*((0.033)^2-(0.026)^2);   % in m^2
    g   = 9.81;                         % Gravity in m/s^2
    rho = 1380;                         % Density in kg/m^3           
    
    % Create the space for stiffness and force matrix
    
    K  = zeros(3*length(mesh.x));   
    f  = zeros(3*length(mesh.x),1);
    
    % Specifying the external forces to the given points in the truss figure
    
    f(length(f))    = f(length(f)) - p;
    f(length(f)-3)  = f(length(f)-3) + (p/2);
    f(length(f)-6)  = f(length(f)-6) + (p/2);
       
    for c           = mesh.conn
        xe          = mesh.x(:,c);                       % Nodal position
        dx          = xe(:,2) - xe(:,1);                 % Difference between 2 nodes (node 2 - node 1)
        t           = dx/norm(dx);                       % Dividing dx by norm(dx) (absolute values)
        Re          = [t' 0 0 0;0 0 0 t'];               % Rotation matrix
        Ke          = Re'*E*A/norm(dx)*[1 -1;-1 1]*Re;   % Element stiffness matrix
        sctr(1:3:6) = 3*c-2;                             % x- degree of freedom is 3*c-2
        sctr(2:3:6) = 3*c-1;                             % y- degree of freedom is 3*c-1
        sctr(3:3:6) = 3*c;                               % z- degree of freedom is 3*c
        w           = rho*g*norm(dx)*A;                  % weight = density*g*length of element*Area
        f(3*c-1)    = f(3*c-1) - (1/2)*w;                   
        K(sctr,sctr)= K(sctr,sctr) + Ke;                 % K matrix (Global)
    end
    
    % Defining the boundary conditions
    
    fixed_dof               = [1*3-2,1*3-1,1*3,2*3-2,2*3-1,2*3,3*3-2,3*3-1,3*3]; 
    K(fixed_dof,:)          = 0;
    K(fixed_dof,fixed_dof)  = eye(length(fixed_dof));
    f(fixed_dof)            = 0;
    
    d = K\f            % Displacement calculations
    sigma       = [];   % Stress (sigma)
    counter     = 1;    % Iteration calculator
    
    for c                   = mesh.conn
        xe                  = mesh.x(:,c);          % Nodal position
        dx                  = xe(:,2) - xe(:,1);    % Difference between 2 nodes (node 2 - node 1)
        t                   = dx/norm(dx);          % Dividing dx by norm(dx) (absolute values)
        Re                  = [t' 0 0 0;0 0 0 t'];  % Rotation matrix
        sctr(1:3:6)         = 3*c-2;                             % x- degree of freedom is 3*c-2
        sctr(2:3:6)         = 3*c-1;                             % y- degree of freedom is 3*c-1
        sctr(3:3:6)         = 3*c;                               % z- degree of freedom is 3*c
        
        original_length     = norm(dx);                                     % Original length
        de                  = d(sctr);                                      % Gather - element displacement
        intermediate_step   = Re*de;         
        delta_L             = intermediate_step(2)-intermediate_step(1);    % Change in length (delta)
        epsilon             = delta_L/original_length;                      % Strain
        
        sigma(counter)      = E*epsilon;
        counter             = counter + 1;
        
    end
    
end

