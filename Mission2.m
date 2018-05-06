clear;
load gauntlet_map/map.mat
%%
clf;
hold on;
while numel(data) > 50
    [x1,x2,y1,y2,inliers] = ransacLine(data,200,0.5,3);
    plot([x1 x2],[y1 y2],'r');
    plot(inliers(:,1),inliers(:,2),'go');
    data = setdiff(data,inliers,'rows');
    numel(data)
end
plot(data(:,1),data(:,2),'b*');
hold off;
%%

%%
function [x1,x2,y1,y2,inliers] = ransacLine(data,num_runs,d_perpendicular, d_gap)
    x = data(:,1);
    y = data(:,2);
    combinations = datasample(nchoosek(1:length(x),2),num_runs);
    calculated_lines = table();
    for combination = combinations'
        rand_1 = combination(1); rand_2 = combination(2);
        x1=x(rand_1);y1=y(rand_1);x2=x(rand_2);y2=y(rand_2);
        points = [x1 y1 x2 y2];
        inliers = 0;
        inlier_points = [];
        for i = 1:size(x,1)
            tx = data(i,1);
            ty = data(i,2);
            distance = abs((y2-y1).*tx - (x2-x1).*ty + x2*y1 - y2*x1)./ ...
                sqrt((y2-y1)^2 + (x2-x1)^2);
            if distance < d_perpendicular
                inliers = inliers + 1;
                inlier_points = [inlier_points;data(i,:)];
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
        current_line = [];
        for i_d = 1:length(dists)
            di = dists(i_d);
            if di < d_gap
                current_line = [current_line ; s_inliers(i_d, :)];
            else
                if ~isempty(current_line)
                    clusters{end+1} = current_line;
                end
                current_line = [];
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