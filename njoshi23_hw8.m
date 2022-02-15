function[N, dNdp, dNdx] = njoshi23_hw8(xe, p)

nne = length(xe);

if nne == 3
    [N, dNdp] = shape_tri3(p);
elseif nne == 6
    [N, dNdp] = shape_tri6(p);
end

J = xe*dNdp;
dNdx = dNdp/J;

end
    
function[N, dNdp] = shape_tri3(p)

N = [p(1);p(2);1-p(1)-p(2)];
dNdp =[1, 0; 0, 1; -1, -1];

end

function[N, dNdp] = shape_tri6(p)

N = [p(1)*(2*p(1)-1);
     p(2)*(2*p(2)-1);
     (1-p(1)-p(2))*(2*(1-p(1)-p(2))-1);
     4*p(1)*p(2);
     4*p(2)*(1-p(1)-p(2));
     4*p(1)*(1-p(1)-p(2))];

dNdp = [4*p(1) - 1,0;
        0,4*p(2) - 1;
        4*p(1)-3+4*p(2),-3+4*p(1)+4*p(2);
        4*p(2), 4*p(1);
        -4*p(2), 4-4*p(1)-8*p(2);
        4-8*p(1)-4*p(2), -4*p(1)];     

end