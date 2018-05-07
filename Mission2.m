clear;
clear global
%%
load gauntlet_map/map.mat
data = data .* 0.0254;
%%
[inliers_circle, inliers_lines] = org_points(data,50,0.00127,0.0762);

figure(fig1)
clf
range_min = -0.762;
range_max = 2.032;
step = 0.05;
[X,Y] = meshgrid([range_min:step:range_max],[range_min:step:range_max]);
Z = generate_scalar_field(inliers_circle,inliers_lines,X,Y,7,1);
Z = reshape(Z,size(X));

hold on
contour(X,Y,Z)

[gx, gy] = gradient(Z);
gx = cap(gx,5);
gy = cap(gy,5);
plot(inliers_circle(:,1),inliers_circle(:,2),'b*');
plot(inliers_lines(:,1),inliers_lines(:,2),'r*');
quiver(X,Y,gx,gy)
hold off

figure(fig2)
clf
surf(X,Y,Z)
title("the real mount doom")
zlabel("degree of terror")
xlabel("wimp zone")
ylabel("hero zone")


% clf;
% hold on;
% plot(data(:,1),data(:,2),'b*');
% plot(inliers_circle(:,1),inliers_circle(:,2),'g*');
% plot(inliers_lines(:,1),inliers_lines(:,2),'r*');
% hold off;
%%
function [x1,x2,y1,y2,inliers] = ransacLine(data,num_runs,d_perpendicular, d_gap)
    if size(data,1) == 1
        x = data(1);
        y = data(2);
        x1 = x; x2 = x; y1 = y; y2 = y;
        inliers = [x y];
        return
    end
    combinations = nchoosek(1:size(data,1),2);
    combinations = datasample(combinations,min(num_runs,size(combinations,1)),1);
%     combinations = nchoosek(1:length(x),2);
    calculated_lines = table();
    for combination = combinations'
        p1 = data(combination(1),:);
        p2 = data(combination(2),:);
        
        t_hat = (p2 - p1)/norm(p2 - p1);
        n_hat = cross([t_hat 0], [0 0 1]);
        n_hat(3) = [];
        
        inliers = 0;
        inlier_points = [];
        for k = 1:size(data,1)
            ptest = data(k,:);
            d = ptest - p1;
            b = abs(dot(n_hat,d));
            if b < d_perpendicular
                inliers = inliers + 1;
                inlier_points = [inlier_points;ptest];
            end
        end
        % create line to find gaps
        clear x2 y2
        [x1, ~, y1, ~] = calc_endpoints(inlier_points);
        % find gaps
        dist = [];
        for inlier = inlier_points'
            tx = inlier(1);
            ty = inlier(2);
            dist = [dist; sqrt((tx-x1)^2 + (ty-y1)^2)];
        end
        [dist, ii] = sort(dist);
        s_inliers = inlier_points(ii,:);
        dists = [0 diff(dist)'];
        % cluster points
        clusters = {};
        current_cluster = [];
        for i_d = 1:length(dists)
            di = dists(i_d);
            if di < d_gap
                current_cluster = [current_cluster ; s_inliers(i_d, :)];
            else
                if ~isempty(current_cluster)
                    clusters{end+1} = current_cluster;
                end
                current_cluster = [];
            end
            if i_d == length(dists)
                clusters{end+1} = current_cluster;
            end
        end
        
        % choose cluster with most points
        clusters_n = [];
        for i = 1:length(clusters)
            clusters_n(end+1) = size(clusters{i},1);
        end
        
        if isempty(clusters_n)
            % skip appending line if there are no clusters
            continue
        end
        
        % recalculate line from largest cluster
        [inliers,clusters_n_max] = max(clusters_n);
        inlier_points = clusters{clusters_n_max};
        [x1, x2, y1, y2] = calc_endpoints(inlier_points);
        points = [x1, y1, x2, y2];
        
        % return variable to access data outside loop
        stuff = {points inliers {inlier_points}};
        calculated_lines = [calculated_lines; stuff];
    end
    [~,I] = max(calculated_lines.b2);
    x1 = calculated_lines.b1(I,1);
    y1 = calculated_lines.b1(I,2);
    x2 = calculated_lines.b1(I,3);
    y2 = calculated_lines.b1(I,4);
    inliers = cell2mat(calculated_lines.b3(I));
end

% function [x1,x2,y1,y2,inliers] = ransacLine(data,num_runs,d_perpendicular, d_gap)
function [center, radius, inliers] = ransacCircle(data, num_runs, d)
    x = data(:,1);
    y = data(:,2);
    combinations = datasample(nchoosek(1:length(x),3),num_runs);
    calculated_circles = table();
    for combination = combinations'
        rand_1 = combination(1); rand_2 = combination(2); rand_3 = combination(3);
        x1=x(rand_1);y1=y(rand_1);x2=x(rand_2);y2=y(rand_2);x3=x(rand_3);y3=y(rand_3);

        A = [...
            (2*x2)-(2*x1) (2*y2)-(2*y1);...
            (2*x3)-(2*x2) (2*y3)-(2*y2)...
        ]; 

        B=[...
            x2^2 - x1^2 + y2^2 - y1^2;...
            x3^2 - x2^2 + y3^2 - y2^2 ...
        ];
        center = linsolve(A,B);
        r = sqrt((x1-center(1,:)).^2 + (y1-center(2,:)).^2);
        inliers = 0;
        inlier_points = [];
        for j = 1:size(x,1)
            xpoint=x(j); ypoint=y(j);
            rpoint = sqrt((xpoint-center(1,:)).^2 + (ypoint-center(2,:)).^2);
            if rpoint < r+d && rpoint > r-d
                inlier_points = [inlier_points; xpoint ypoint];
                inliers = inliers + 1;
            end
        end
        stuff = {center' r inliers {inlier_points}};
        calculated_circles = [calculated_circles; stuff];
    end
    [~,I] = max(calculated_circles{:,3});
    center = calculated_circles.b1(I,1:2);
    radius = calculated_circles.b2(I);
    inliers = calculated_circles.b4{I};
end

function [x1, x2, y1, y2] = calc_endpoints(points) 
    % create line to find gaps
    [y1, i1] = min(points(:,2));
    [y2, i2] = max(points(:,2));
    if y1 ~= y2
       x1 = points(i1,1);
       x2 = points(i2,1);
    else
        [x1, i1] = min(points(:,1));
        [x2, i2] = max(points(:,1));
        y1 = points(i1,2);
        y2 = points(i2,2);
    end
end

function [inliers_circle, inliers_lines] = org_points(data,num_runs,d_perpendicular, d_gap)
    data = data;
    found_circle = false;
    l_in = []; c_in = [];
    while numel(data) > 0
        [x1,x2,y1,y2,inliers_line] = ransacLine(data,num_runs,d_perpendicular, d_gap);

        if ~found_circle
            [center,radius,inliers_circle] = ransacCircle(data, 500, d_perpendicular);
        end

        if found_circle || size(inliers_line,1) > size(inliers_circle,1)
            inliers = inliers_line;
            l_in = [l_in; inliers_line];
        else
            found_circle = true;
            inliers = inliers_circle;
            c_in = inliers_circle;
        end
        data = setdiff(data,inliers,'rows');
    end
    inliers_circle = c_in;
    inliers_lines = l_in;
end

function out = cap(x, threshold)
    x(x>0 & x>threshold) = threshold;
    x(x<0 & abs(x)>threshold) = -threshold;
    out = x;
end

function z = generate_scalar_field(attract_points, repell_points, grid_x, grid_y, attract_power, repell_power)
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
            z_i = z_i - attract_power*log(sqrt((x-x_i)^2 + (y-y_i)^2));
        end
        
        for i_repell = 1:N_repell
            repell_point = repell_points(i_repell,:);
            x_i = repell_point(1); y_i = repell_point(2); 
            z_i = z_i + repell_power*log(sqrt((x-x_i)^2 + (y-y_i)^2));
        end
        
        z(i_grid) = z_i;
    end
end