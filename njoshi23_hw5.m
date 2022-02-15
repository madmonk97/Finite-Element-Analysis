function [N, dNdp] = njoshi23_hw5(p, nne)
if nne == 2
    [N,dNdp] = shape2(p);
elseif nne == 3
    [N,dNdp] = shape3(p);
elseif nne == 4
    [N,dNdp] = shape4(p);
end
end

function [N,dNdp] = shape2(p)
N   = (1/2)*[1-p;1+p];
dNdp= [-0.5;0.5];
end

function [N,dNdp] = shape3(p)
N   = [(-1/2)*p*(1-p) ; (1-p)*(1+p) ; (1/2)*p*(1+p)];
dNdp= [p-0.5 ; -2*p ; p+0.5];
end

function [N,dNdp] = shape4(p)
N   = [(-9/16)*(p+1/3)*(p-1/3)*(p-1) ; (27/16)*(p+1)*(p-1/3)*(p-1) ; (-27/16)*(p+1)*(p+1/3)*(p-1) ; (9/16)*(p+1)*(p+1/3)*(p-1/3)];
dNdp= [(-9/16)*(3*p^2-2*p-1/9) ; (27/16)*(3*p^2-2/3*p-1) ; (-27/16)*(3*p^2+2/3*p-1) ; (9/16)*(3*p^2+2*p-1/9)];
end


