function [f] = code_7(mesh, b)
% Defines default inputs for testing.
if nargin == 0
    ne = 4;
    nn = 2*ne + 1;
    mesh.x = linspace(-2, 2, nn);
    mesh.conn = [1:2:nn-2; 2:2:nn-1; 3:2:nn];
    b = @(x) exp(-x.^2);
end
% Compute nodal force vector.
qpts = quadrature(2);
f = zeros(length(mesh.x),1);

for c = mesh.conn
    xe = mesh.x(:,c);
    b_tot = 0;

    for q = qpts
        [N,dNdp] = shape3(q(1));
        J = xe*dNdp;
        x = xe*N;
        b_tot = b_tot + b(x)*N*q(2)*J;
    end
    f(c) = f(c)+b_tot;
end

end


function [qpts] = quadrature(n)
%   QUADRATURE
%     quadrature(n) returns a quadrature table for a rule with n
%     integration points.  The first row of the table gives the quadrature
%     point location, and the second gives the quadrature weights.

u = 1:n-1;
u = u./sqrt(4*u.^2 - 1);

A = zeros(n);
A(2:n+1:n(n-1)) = u;
A(n+1:n+1:n^2-1) = u;

[v, x] = eig(A);
[x, k] = sort(diag(x));
qpts = [x'; 2*v(1,k).^2];
end

function [N, dNdp] = shape3(p)
 N = [(-1/2)*p*(1-p) ; (1-p)*(1+p) ; (1/2)*p*(1+p)];
 dNdp = [p - 0.5 ; -2*p ; p + 0.5];
end