function [e1, e2] = code_6(mesh, d)

clc; clear;

% Defines default inputs for testing.
    
    if nargin == 0
        ne = 4;
        nn = 2*ne + 1;
        mesh.x = linspace(-2, 2, nn);
        mesh.conn = [1:2:nn-2; 2:2:nn-1; 3:2:nn];
        d = (sin(2*pi*mesh.x.^2).*exp(-mesh.x.^2))';
        E = 10; % reference
        A = 10; % reference
    end
    
% Compute e1 and e2
quadr_tbl= [-sqrt(3/7+2/7*sqrt(6/5)), -sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5)); (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];
K        = zeros(length(mesh.x)); 

for c   = mesh.conn
    xe  = mesh.x(:,c);
    Ke  = zeros(length(c));
    
    for q = quadr_tbl
        [~,dNdp] = shape3(q(1));    % q(1)
        J = xe*dNdp;                % dNdx ??
        B = dNdp/J;
        Ke = Ke + B*E*A*B'*J*q(2);  % q(2) is weight
    end
    
    sctr = c;
    K(sctr,sctr) = K(sctr,sctr) + Ke;

end

    function [N,dNdp] = shape3(p)
        N    = [-0.5*p*(1-p) ; (1-p)*(1+p) ; 0.5*p*(1+p)];
        dNdp = [p-0.5 ; -2*p ; p+0.5];
    end

end