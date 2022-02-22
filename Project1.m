function d = Project1(ne,eltype)
%% Project 1

%Defining key values
dia = 0.05; %meters
E = 70*10^9; %Pascals
rho = 2700; %kg/m^3
m = 1; %kg
r = 1; %m
A = (0.05/2)^2*pi; %m^2
omega = 2000 * 2 * pi / 60; %rad/s

% Centripetal Force is m*v^2/r; v is r*omega
v = r*omega; %m/s
Fc = m*v^2/r; %N


%All of this is set up for linear elements
if nargin == 0
    ne = 2; %number of elements is 2
    eltype = 1; %linear
end
nn = ne + 1; %calculated number of nodes
mesh.x = linspace(0, 1, nn); %creating mesh
mesh.conn = [1:1:nn-1 ; 2:1:nn]; %connectivity matrix

%Preallocation
K = zeros(length(mesh.x));
f = zeros(length(mesh.x),1);

%% Stiffness Matrix Creation

if eltype == 1
   quadr_tbl = [-1/sqrt(3), 1/sqrt(3); 1, 1];
elseif eltype == 2
   quadr_tbl = [-sqrt(3/5), 0, sqrt(3/5); 5/9, 8/9, 5/9];
elseif eltype == 3
   quadr_tbl = [-sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5));(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];             
end
for c = mesh.conn
   
    xe = mesh.x(:,c);
    Ke = zeros(length(c));
    
    for q = quadr_tbl
        if eltype == 1
            [N,dNdp] = shape2(q(1));
        elseif eltype == 2
            [N,dNdp] = shape3(q(1));
        elseif eltype == 3
            [N,dNdp] = shape4(q(1));
        end
        
        
        J = xe * dNdp;
        dNdx = dNdp / J; %This is B in notes
        
        Ke = Ke + dNdx*E*A*dNdx'*J*q(2);
               
    end
    
    K(c,c) = K(c,c) + Ke; %1D lets us do this
    
end

%% Handling Forces and BC

%Body force loop
for i = 1:length(mesh.x)
   x_i = mesh.x(i);
   b = bodyforce(x_i);
   f(i) = f(i) + b;
end

f(length(f)) = f(length(f)) + m*(r*omega)^2; %Apply the mass to the end

fixdof = 1;
K(fixdof,:) = 0;
K(fixdof,fixdof) = eye(length(fixdof));
f(fixdof) = 0;

d = K\f;

%% Strain calcs

xx = [];
strain = [];
uu = [];
for c = mesh.conn
    
    de = d(c)';
    xe = mesh.x(:,c);
    for p = linspace(-1,1,20)
       if eltype == 1
          [N, dNdp] = shape2(p);
       elseif eltype == 2
          [N, dNdp] = shape3(p);
       elseif eltype == 3
          [N, dNdp] = shape4(p);
       end
       
       J = xe*dNdp;
       dNdx = dNdp/J;
       
       xx(length(xx)+1) = xe*N;
       uu(length(uu)+1) = de*N;
       strain(length(strain) + 1) = de*dNdx;
       
    end
    
end

end


%% Extra needed functions
function [N, dNdp] = shape2(p)
 N = (1/2)*[1-p;1+p];
 dNdp = [-0.5;0.5];
end
function [N, dNdp] = shape3(p)
 N = [(-1/2)*p*(1-p) ; (1-p)*(1+p) ; (1/2)*p*(1+p)];
 dNdp = [p - 0.5 ; -2*p ; p + 0.5];
end
function [N, dNdp] = shape4(p)
 N = [(-9/16)*(p+1/3)*(p-1/3)*(p-1) ; (27/16)*(p+1)*(p-1/3)*(p-1) ; (-27/16)*(p+1)*(p+1/3)*(p-1) ; (9/16)*(p+1)*(p+1/3)*(p-1/3)];
 dNdp = [(-9/16)*(3*p^2-2*p-1/9) ; (27/16)*(3*p^2-2/3*p-1) ; (-27/16)*(3*p^2+2/3*p-1) ; (9/16)*(3*p^2+2*p-1/9)];
end

function b = bodyforce(x)
    rho = 2700; %kg/m^3
    dia = 0.05; %meters
    A = pi*(dia/2)^2; %m^2
    omega = 2000 * 2 * pi /60; %rad/s
    v = omega.*x; %m/s
    b = rho*A.*v.^2./x;
end