clear all

%%
x_neato = -2;
y_neato = 0;

x_bob = 0.33274;
y_bob = 1.8288;

syms x y
f = -log(sqrt((x -x_bob)^2 + (y-y_bob)^2));
g = gradient(f,[x,y]);

g_neato = subs(g,[x,y],{x_neato,y_neato});
double(norm(g_neato))

%%
sub_bump = rossubscriber('/bump');
sub_encoder = rossubscriber('/encoders');
v_pub = rospublisher('/raw_vel');
v_msg = rosmessage(v_pub);
d = 0.234;
v_max = 0.15;

data_encoder = [];
data_theta = [0];
data_r = [0 0];

figure(r_fig);
clf
hold on
ax = gca;
plot(x_bob, y_bob, "r*");
q_grad = quiver(0,0,0,0,"r")
q_heading = quiver(0,0,0,0,"b")
axis([-1 3 -1 3])
title("Finding Bob: Mission 1")

"starting..."
bump = 0;
while ~bump
    bumpd = receive(sub_bump);
    bump = any(bumpd.Data);
    
    pause(0.1)
    
    encoder = receive(sub_encoder);
    encoder = encoder.Data;
    
    % skip first + second cycle
    if size(data_encoder,1) == 0
        data_encoder(1,:) = encoder';
        continue
    elseif size(data_encoder,1) == 1
        data_encoder(2,:) = encoder';
        continue
    end
    
    data_encoder = [data_encoder(2,:) ; encoder'];
    v_wheels = diff(data_encoder);
    v = mean(v_wheels);
    
    w = rad2deg((v_wheels(2) - v_wheels(1))/d);
    theta = data_theta(end) + w;
    T_hat = [cosd(theta) sind(theta)];
    drdt = v*T_hat;
    r = data_r(end,:) + drdt;
    
    data_theta(end+1) = theta;
    data_r(end+1,:) = r;
    
    plot(ax, r(1),r(2),"b*");
    
    g_neato = subs(g,[x,y],{r(1),r(2)});
    g_neato = double(g_neato);
    
    u = [g_neato'];
    v = [T_hat];
    phi = -(atan2d(u(1)*v(2)-u(2)*v(1),u(1)*v(1)+u(2)*v(2)));
    
    a = (1/180)*phi;
    b = (-1/180)*phi + 1;
    rotation = ((180-phi)/abs(180-phi))*v_max*a;
    forward = v_max*b;
    
    vr = rotation + forward;
    vl = -rotation + forward;
    sendVel(vl,vr,v_pub,v_msg);
    
    set(q_grad,...
        'xdata',r(1),...
        'ydata',r(2),...
        'udata',r(1)+g_neato(1),...
        'vdata',r(2)+g_neato(2))
    
    set(q_heading,...
        'xdata',r(1),...
        'ydata',r(2),...
        'udata',r(1)+T_hat(1),...
        'vdata',r(2)+T_hat(2))
end
hold off

sendVel(0,0,v_pub,v_msg);
% tftree = rostf;
% odom_pose = getTransform(tftree, 'odom', 'base_link')
% trans = odom_pose.Transform.Translation;
% quat = odom_pose.Transform.Rotation;
% rot = quat2eul([quat.W quat.X quat.Y quat.Z])
%%

function sendVel(vl, vr, pub, msg)
    sprintf("VL: %.4d, VR: %.4d", vl, vr);
    msg.Data = [vl vr];
    send(pub, msg);
end