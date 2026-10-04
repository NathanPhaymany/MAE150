clear all;
close all;
clc;
format long;

name = 'Nathan Phaymany';
id = '20016227';
hw_num = 1;

%----------------------------------------------------------------------------------------------------------
% Basic Functions
plot_points2 = @(ax, list_of_points, cs) plot(ax, list_of_points(1, :), list_of_points(2, :), [cs, '-']);
plot_points3 = @(ax, list_of_points, cs) plot3(ax, list_of_points(1, :), list_of_points(2, :), list_of_points(3, :), [cs, '-']);
e = exp(1);


%----------------------------------------------------------------------------------------------------------
% Part 1

S = @(sx, sy) [sx 0 0; 0 sy 0; 0 0 1];
T = @(tx, ty) [eye(3, 2) [tx ty 1]'];
R = @(theta) [cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; 0 0 1];

p = ones(3, 361);
P1 = ones(3, 361);
P2 = ones(3, 361);
P3 = ones(3, 361);

for s = 0:360
    x = cosd(s);
    y = sind(s);
    p(:, s+1) = [x; y; 1];
    
    P1(:, s+1) = R(-45)*S(2,1)*p(:, s+1);
    P2(:, s+1) = R(-45)*S(4,1)*p(:, s+1);
    P3(:, s+1) = R(-45)*S(6,1)*p(:, s+1);
end

f1 = figure(1);
ax1 = axes('Parent', f1);
hold(ax1, 'on')

plot_points2(ax1, p, 'k')
plot_points2(ax1, P1, 'm')
plot_points2(ax1, P2, 'b')
plot_points2(ax1, P3, 'r')
xlim(ax1, [-5 5])
ylim(ax1, [-5 5])
title(ax1, "2D Transformation of a parametric curve")
xlabel(ax1, 'x')
ylabel(ax1, 'y')
legend(ax1, 'circle radius = 1', 'major axis radius = 2', 'major axis radius = 4', 'major axis radius = 6')
grid(ax1, 'on')
p1 = 'see figure 1';


%----------------------------------------------------------------------------------------------------------
% Part 2
clear S T R p x y;
S = @(sx, sy, sz) [sx,  0,  0,   0; ...
                    0, sy,  0,   0; ...
                    0,  0, sz,   0; ...
                    0,  0,  0,   1];     
           
% Translation matrix        
T = @(tx, ty, tz) [eye(4,3) [tx, ty, tz, 1]'];

% Rotation matrix around z axis: positive thetaz in degree counterclockwise
Rz = @(thetaz) [cosd(thetaz),  -sind(thetaz),   0,   0; ...
                sind(thetaz),   cosd(thetaz),   0,   0; ...
                            0,             0,   1,   0; ...
                            0,             0,   0,   1];

% Rotation matrix around x axis: positive thetax in degree counterclockwise
Rx = @(thetax) [1,              0,             0,    0; ...
                0,   cosd(thetax),  -sind(thetax),   0; ...
                0,   sind(thetax),   cosd(thetax),   0; ...
                0,              0,              0,   1];

% Rotation matrix around y axis: positive thetay in degree counterclockwise
Ry = @(thetay) [cosd(thetay),   0,   sind(thetay),   0; ...
                           0,   1,              0,   0; ...      
               -sind(thetay),   0,   cosd(thetay),   0;...
                           0,   0 ,             0,   1];


s = linspace(0, 12*pi, 1000);
x = sin(s) .* (exp(cos(s)) - 2*cos(4*s) - sin(0.1 * s).^2);
y = -cos(s) .* (exp(cos(s)) - 2*cos(4*s) - sin(0.1 * s).^2);
z = zeros(size(x));

p = [x; y; z; ones(size(x))];


f2 = figure(2);
ax2 = axes('Parent', f2);

hold(ax2, 'on')
h2 = plot_points3(ax2, nan(4, 1), 'b');
grid(ax2, 'on')
xlabel(ax2, 'x (cm)')
ylabel(ax2, 'y (cm)')
zlabel(ax2, 'z (cm)')
xlim(ax2, [-10 10])
ylim(ax2, [-15 0])
zlim(ax2, [-5 5])
axis(ax2, 'manual')
view(ax2, 3)
titleHandle = title(ax2, '3D CAD transformation: time = 0.0 s');

yoff = 0;

for t = 0:0.1:10 
    transform_points = T(0, -yoff, 0)*Ry(30*sin(3*pi*t/2))*p;

    set(h2, 'XData', transform_points(1, :), ...
            'YData', transform_points(2, :), ...
            'ZData', transform_points(3, :));
    titleHandle.String = sprintf('3D CAD transformation: time = %.1f s', t);
    yoff = yoff + 0.1;
    pause(0.1);
end

