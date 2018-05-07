
%Change these ranges and point spacings to make nice figures
range_min = -4;
range_max = 4;
[X,Y] = meshgrid([range_min:0.1:range_max],[range_min:0.1:range_max]);
[X1,Y1] = meshgrid([[range_min:0.5:-0.3] [0.3:0.5:range_max]],[[range_min:1:-1] [1:1:range_max]]);

%Define your function below.  Don't forget to specify element by element
%operations.
syms x y
xi = 1
yi = 2
% f = 12-(x-2).^2 -(y-4).^2;
% f = sin(x) + cos(y);
f = -log(sqrt((x -xi)^2 + (y-yi)^2))

%This takes the gradient of your function
g = gradient(f,[x,y])

%These lines of code make nice plots.  Look at documentation to adjust
%plotting parameters as you like.
hold off
contour(X,Y,subs(f,[x,y],{X,Y}))

G1 = subs(g(1),[x,y],{X1,Y1});
G2 = subs(g(2),[x,y],{X1,Y1});
hold on
quiver(X1,Y1,G1,G2)

%%
clf
range_min = -4;
range_max = 4;
[X,Y] = meshgrid([range_min:0.1:range_max],[range_min:0.1:range_max]);
[X1,Y1] = meshgrid([[range_min:0.5:-0.3] [0.3:0.5:range_max]],[[range_min:1:-1] [1:1:range_max]]);
Z = generate_scalar_field([0 0 ; -1 -1],[2 2; 3 3],X,Y);
Z = reshape(Z,size(X));

hold on
contour(X,Y,Z)

[gx, gy] = gradient(Z);

gx = cap(gx,1);
gy = cap(gy,1);
quiver(X,Y,gx,gy)
hold off

%%
function out = cap(x, threshold)
    x(x>0 & x>threshold) = threshold;
    x(x<0 & abs(x)>threshold) = -threshold;
    out = x;
end

function z = generate_scalar_field(attract_points, repell_points, grid_x, grid_y)
    N_grid = numel(grid_x);
    N_attract = size(attract_points,1);
    N_repell = size(repell_points,1);
    z = zeros(N_grid,1);
    for i_grid = 1:N_grid
        x = grid_x(i_grid); y = grid_y(i_grid);

        z_i = 0;
        for i_attract = 1:N_attract
            attract_point = attract_points(i_attract,:);
            x_i = attract_point(1); y_i = attract_point(2); 
            z_i = z_i - log(sqrt((x-x_i)^2 + (y-y_i)^2));
        end
        
        for i_repell = 1:N_repell
            repell_point = repell_points(i_repell,:);
            x_i = repell_point(1); y_i = repell_point(2); 
            z_i = z_i + log(sqrt((x-x_i)^2 + (y-y_i)^2));
        end
        
        z(i_grid) = z_i;
    end
end