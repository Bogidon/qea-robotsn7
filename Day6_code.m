clear all;
load scan4.mat;
%%
theta_clean = theta(find(r));
r_clean = r(find(r));
[x,y] = pol2cart(deg2rad(theta_clean), r_clean);
data = [x,y];
%%
num_runs = 200;
d = 0.02;
%%
combinations = datasample(nchoosek(1:length(x),2),num_runs);
calculated_lines = table();
clf;hold on;plot(x,y,'cs');
for combination = combinations'
    rand_1 = combination(1); rand_2 = combination(2);
    x1=x(rand_1);y1=y(rand_1);x2=x(rand_2);y2=y(rand_2);
    points = [x1 y1 x2 y2];
    inliers = 0;
    inlier_points = [];
    for i = 1:size(x,1)
        distance = abs((points(:,4)-points(:,2)).*data(i,1) - (points(:,3)-points(:,1)).*data(i,2) + points(:,3)*points(:,2) - points(:,4)*points(:,1))./ ...
            sqrt((points(:,4)-points(:,2))^2 + (points(:,3)-points(:,1))^2);
        if distance < d
            inliers = inliers + 1;
            inlier_points = [inlier_points;data(i,:)];
        end
    end
    calculated_lines = [calculated_lines;{points inliers inlier_points}];
end
[~,I] = max(calculated_lines.b2);
px1 = calculated_lines.b1(I,1);
py1 = calculated_lines.b1(I,2);
px2 = calculated_lines.b1(I,3);
py2 = calculated_lines.b1(I,4);
ins = cell2mat(calculated_lines.b3(I));
% plot(ins(:,1),ins(:,2),'r*')
miny = min(ins(:,2));
maxy = max(ins(:,2));
slope = (py2 - py1)/(px2 -px1);
b = py1 - (slope * px1);syms t;
xmin = (miny-b)/slope;xmax = (maxy-b)/slope;
fplot(slope*(t-px1)+py1,[xmin,xmax],'-r');
% plot(calculated_lines(I,1:2:3),calculated_lines(I,2:2:4),'r');
hold off;