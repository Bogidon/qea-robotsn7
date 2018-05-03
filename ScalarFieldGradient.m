
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