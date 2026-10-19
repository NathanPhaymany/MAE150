clear all; 
close all;
clc;
format long;

name = 'Nathan Phaymany';
id = '20016227';
hw_num = 3;

%% Problem 1

syms t
N1 = 2*t^3 - 3*t^2 + 1;
N2 = t^3 - 2*t^2 + t;
N3 = -2*t^3 + 3*t^2;
N4 = -t^2 + t^3;

x = [-5 0; 0 4; 4 7; 7 9];
y = [0 5; 5 3; 3 -2; -2 1];
dydx = [1 -2; -2 1; 1 -5; -5 0];
D = [-0.4 0.4];

figure(1); hold on;
cs = 'krbgcmy';

for i = 1:4
    % Initial Curves
    C{1, i} = (x(i, 2) - x(i, 1)) * t + x(i, 1);
    C{2, i} = y(i, 1) * N1 ...
        + dydx(i, 1) * (x(i, 2) - x(i, 1)) * N2 ...
        + y(i, 2) * N3 ...
        + dydx(i, 2) * (x(i, 2) - x(i, 1)) * N4;
    
    xc(i, :) = subs(C{1, i}, t, 0:0.1:1);
    yc(i, :) = subs(C{2, i}, t, 0:0.1:1);

    plot(xc(i, :), yc(i, :), cs(i));
    lg{i} = sprintf('Curve #%i', i);

    % Offset curves
    s = sqrt(diff(C{1, i}, t, 1)^2 + diff(C{2, i} , t, 1)^2);
    % x-component of unit tangent
    T{1, i} = 1/s * diff(C{1, i}, t, 1);
    % y-component of unit tanget
    T{2, i} = 1/s * diff(C{2, i}, t, 1);
    % x-component of unit normal
    N{1, i} = 1/s * -diff(C{2, i}, t, 1);
    % y-component of unit normal
    N{2, i} = 1/s *  diff(C{1, i}, t, 1);

    for n = 1:length(D)
        % x-component of offset
        O{1, n} = C{1, i} + D(n) * N{1, i};
        % y-component of offset
        O{2, n} = C{2, i} + D(n) * N{2, i};
        % evaluate x- and y-values of the offset for given range of t
        x_off(n, :) = subs(O{1, n}, t, 0:0.01:1);
        y_off(n, :) = subs(O{2, n}, t, 0:0.01:1);

        plot(x_off(n, :), y_off(n, :), cs(4+n));
        lg{4+n} = sprintf('D = %.1f',D(n));
    end
end

legend(lg);
xlabel('x');
ylabel('y');
title('Problem 1');
box on; grid on;

p1a = 'See figure 1'
p1b = 'The offset curve when D = -0.4 is smooth. However, the offset curve when D = 0.4 is not smooth as the curve experienced undercutting (around x = 7.5) since the offset distance is too large to smoothly round the curve.'

%% Problem 2

clear t;
syms theta beta L

dwell = L;
harmonic = L/2 * (1-cos(pi*theta/beta));
cycloidal = L * (theta/beta - sin(2*pi*theta/beta)/(2*pi));
poly345 = L * (10/beta^3 * theta^3 - 15/beta^4 * theta^4 + 6/beta^5 * theta^5);

r_follower = 0.006;
r_prime = 0.02;
omega = 50*2*pi/60;
r_base = r_prime - r_follower;
dt = deg2rad(2);
e = 0;
R0 = sqrt(r_prime^2-e^2);

% Section 1: 
n = 1;
cam(n).c_type = poly345;        % curve type
cam(n).L_range = [0, 10];          % [start_value end_value]          
cam(n).theta_range = deg2rad([0, 120]);     % [start_theta end_theta]

% Section 2:  
n = 2;
cam(n).c_type = dwell;
cam(n).L_range = [10, 10];
cam(n).theta_range = deg2rad([120, 130]);

% Section 3:  
n = 3;
cam(n).c_type = cycloidal;
cam(n).L_range = [10, 20];
cam(n).theta_range = deg2rad([130, 210]);

% Section 4:  
n = 4;
cam(n).c_type = dwell;
cam(n).L_range = [20, 20];
cam(n).theta_range = deg2rad([210, 230]);

% Section 5:  
n = 5;
cam(n).c_type = harmonic;
cam(n).L_range = [20, 0];
cam(n).theta_range = deg2rad([230, 330]);

% Section 6:  
n = 6;
cam(n).c_type = dwell;
cam(n).L_range = [0, 0];
cam(n).theta_range = deg2rad([330, 360]);

for n = 1:length(cam) 

    % Set up theta with radial resolution dt
    if n == length(cam)
        t = cam(n).theta_range(1) : dt : cam(n).theta_range(2); 
    else
        t = cam(n).theta_range(1) : dt : cam(n).theta_range(2)-dt; 
    end

    y = cam(n).L_range(1) ...
        + subs(cam(n).c_type, {theta, beta, L}, ...
        {t - t(1), diff(cam(n).theta_range), diff(cam(n).L_range)});
        
    % Compute velocity function v = dy/dtime = dy/dtheta * omega
    v = subs(omega * diff(cam(n).c_type, theta, 1), {theta, beta, L}, ...        
        {t - t(1), diff(cam(n).theta_range), diff(cam(n).L_range)});

    % Compute acceleration function  a = d2y/dtime2 = d2y/dtheta2 * omega^2
    a = subs(omega^2 * diff(cam(n).c_type, theta, 2), {theta, beta, L}, ...        
        {t - t(1), diff(cam(n).theta_range), diff(cam(n).L_range)});
    
    % Compute jerk function j = d3y/dtime3 = d3y/dtheta3 * omega^3
    j = subs(omega^3 * diff(cam(n).c_type, theta, 3), {theta, beta, L}, ...        
        {t - t(1), diff(cam(n).theta_range), diff(cam(n).L_range)});
    
    % Compute pressure angle function with eccentricity e = 0;
    phi = atan((v/omega - e) ./ (y + R0));

    % Evaluate SVAJ and phi at theta values:
    cam(n).theta = t;
    cam(n).dis = double(y);
    cam(n).vel = double(v);
    cam(n).acc = double(a);
    cam(n).jerk = double(j);
    cam(n).phi = double(phi);
end

t = [cam.theta];
dis = [cam.dis];
vel = [cam.vel];
acc = [cam.acc];
jerk = [cam.jerk];
phi = [cam.phi];

% Part a
figure(2);
% Displacement
ax1 = subplot(5, 1, 1);
plot(t, dis);
ylabel({'displacement', '[mm]'});

% Velocity
ax2 = subplot(5, 1, 2);
plot(t, vel);
ylabel({'velocity', '[mm/s]'});

% Acceleration
ax3 = subplot(5, 1, 3);
plot(t, acc);
ylabel({'acceleration', '[mm/s^2]'});

% Jerk
ax4 = subplot(5, 1, 4);
plot(t, jerk);
ylabel({'jerk', '[mm/s^3]'});

% Pressure Angle
ax5 = subplot(5, 1, 5);
plot(t, phi);
ylabel({'\phi', '[degree]'});

ax = [ax1, ax2, ax3, ax4, ax5];
for i = 1:length(ax)
    box(ax(i), 'on');
    grid(ax(i), 'on');
    xlim(ax(i), [0 t(end)]);
end

xlabel('\theta (degree)');

p2a = 'See figure 2'


% Part b
figure(3)

x_pitch = (r_prime + dis) .* sin(t);
y_pitch = (r_prime + dis) .* cos(t);

x_cam = x_pitch - r_follower * sin(t - phi);
y_cam = y_pitch - r_follower * cos(t - phi);

x_base = r_base * sin(t);
y_base = r_base * cos(t);

x_prime = r_prime * sin(t);
y_prime = r_prime * cos(t);

cam_in(1, :) = x_cam;
cam_in(2, :) = y_cam;
pitch_in(1, :) = x_pitch;
pitch_in(2, :) = y_pitch;

cam_in(:, end) = cam_in(:, 1);
pitch_in(:, end) = pitch_in(:, 1);

plot(x_base, y_base, '--g', x_pitch, y_pitch, '-.r', x_cam, y_cam, '-k', x_prime, y_prime, ':g');

hold on; axis equal; box on; grid on;
xlabel('x [mm]');
ylabel('y [mm]');
legend('base', 'pitch', 'cam', 'prime');

p2b = 'See figure 3'
