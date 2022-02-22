skipearly = 1;

if skipearly == 0
%% Tri 3 Elements

dowewantflux = 0;

[d, mesh, xx, qq, p, p2, ~] = Project2(1, 'tri3', 1);
patch(p)
colorbar
axis square
xlabel('x coordinate (mm)')
ylabel('y coordinate (mm)')
title('tri 3 elements, T field, h = 1mm, r = 1mm')


if dowewantflux == 1
    figure
    patch(p2)
    colorbar
    axis square
    xlabel('x coordinate (mm)')
    ylabel('y coordinate (mm)')
    title('tri 3 elements, flux field, h = 1mm, r = 1mm')
    hold on
    quiver(xx(:,1),xx(:,2),qq(:,1),qq(:,2),1,'color','w')
end

%% Quad 4 Elements

figure
[d, mesh, xx, qq, p, p2, ~] = Project2(1, 'quad4', 1);
patch(p)
colorbar
axis square
xlabel('x coordinate (mm)')
ylabel('y coordinate (mm)')
title('quad 4 elements, T field, h = 1mm, r = 1mm')


if dowewantflux == 1
    figure
    patch(p2)
    colorbar
    axis square
    xlabel('x coordinate (mm)')
    ylabel('y coordinate (mm)')
    title('quad 4 elements, flux field, h = 1mm, r = 1mm')
    hold on
    quiver(xx(:,1),xx(:,2),qq(:,1),qq(:,2),1,'color','w')
end

%% Tri 6 Elements

figure
[d, mesh, xx, qq, p, p2, ~] = Project2(1, 'tri6', 1);
patch(p)
colorbar
axis square
xlabel('x coordinate (mm)')
ylabel('y coordinate (mm)')
title('tri 6 elements, T field, h = 1mm, r = 1mm')


if dowewantflux == 1
    figure
    patch(p2)
    colorbar
    axis square
    xlabel('x coordinate (mm)')
    ylabel('y coordinate (mm)')
    title('tri 6 elements, flux field, h = 1mm, r = 1mm')
    hold on
    quiver(xx(:,1),xx(:,2),qq(:,1),qq(:,2),1,'color','w')
end

%% Quad 8 Elements

figure
[d, mesh, xx, qq, p, p2, ~] = Project2(1, 'quad8', 1);
patch(p)
colorbar
axis square
xlabel('x coordinate (mm)')
ylabel('y coordinate (mm)')
title('quad 8 elements, T field, h = 1mm, r = 1mm')


if dowewantflux == 1
    figure
    patch(p2)
    colorbar
    axis square
    xlabel('x coordinate (mm)')
    ylabel('y coordinate (mm)')
    title('quad 8 elements, flux field, h = 1mm, r = 1mm')
    hold on
    quiver(xx(:,1),xx(:,2),qq(:,1),qq(:,2),1,'color','w')
end
end


%% Refinement test
refinementskip = 1;

if refinementskip == 0
sizes = [2, 1,  0.5, 0.25, 0.125];
Ts = zeros(1,length(sizes));
qs = zeros(1,length(sizes));

for i = 1:length(sizes)
    [d, mesh, xx, qq, p, p2, qmax] = Project2(sizes(i), 'quad8', 1);
    Ts(i) = max(d);
    qs(i) = qmax;
end

refinementtbl = [sizes; Ts; qs]
end

%% Second Refinement Test
refinementskip2 = 0;

if refinementskip2 == 0
    
skipr1 = 1;
skipr2 = 1;
if skipr1 == 0
[~, ~, ~, ~, ~, ~, q_max_r1] = Project2LR(0.125, 'quad8', 1);
end
if skipr2 == 0
[~, ~, ~, ~, ~, ~, q_max_r2] = Project2LR(0.125, 'quad8', 2);
end
[~, ~, ~, ~, ~, ~, q_max_r4] = Project2LR(0.125, 'quad8', 4);

max_q_r1 = max(q_max_r1);
max_q_r2 = max(q_max_r2);
max_q_r4 = max(q_max_r4);

e1 = 1;
e2 = 1;
e4 = 1;

result_r1 = [];
result_r2 = [];
result_r4 = [];

sizes = [5 2 1 0.5 0.25 0.125];


if skipr1 == 0 
i = 1;
while abs(e1) > 0.05 && i < 7
    
    [~, ~, ~, ~, ~, ~, q_r1] = Project2(sizes(i), 'quad8', 1);
    e1 = (q_r1 - max_q_r1)/max_q_r1;
    result_r1(end + 1, :) = [q_max_r1, q_r1, e1]';
    i = i + 1;
    
end
end

if skipr2 == 0
i = 1;
while abs(e2) > 0.05 && i < 7
    
    [~, ~, ~, ~, ~, ~, q_r2] = Project2(sizes(i), 'quad8', 2);
    e2 = (q_r2 - max_q_r2)/max_q_r2;
    result_r2(end + 1, :) = [q_max_r2, q_r2, e2]';
    i = i + 1;
    
end
end

i = 1;
while abs(e4) > 0.05 && i < 7
    
    [~, ~, ~, ~, ~, ~, q_r4] = Project2(sizes(i), 'quad8', 4);
    e4 = (q_r4 - max_q_r4)/max_q_r4;
    result_r4(end + 1, :) = [q_max_r4, q_r4, e4]';
    i = i + 1;
    
end
    
end