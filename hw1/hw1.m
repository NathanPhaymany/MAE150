clear all;
close all;
clc;
format long;

name = 'Nathan Phaymany';
id = '20016227';
hw_num = 1;

% Basic Functions
function translation = T(tx, ty)
    transation = [1 0 tx; 0 1 ty; 0 0 1];
end

function rotation = R(theta)
    rotation = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
end

function scaling = S(sx, sy)
    scaling = [sx 0 0; 0 sy 0; 0 0 1];
end


% Part 1

p = ones(3, 361);
P1 = ones(3, 361);
P2 = ones(3, 361);
P3 = ones(3, 361);

for s = 0:360
    x = cos(s);
    y = sin(s);
    p(:, s+1) = [x; y; 1];
    
    P1(:, s+1) = R(-45)*S(2,1)*p(:, s+1);
    P2(:, s+1) = R(-45)*S(4,1)*p(:, s+1);
    P3(:, s+1) = R(-45)*S(6,1)*p(:, s+1);
end

plot(p(:, 1), p(:, 2), '-k', P1(:, 1), P1(:, 2), '-p', P2(:, 1), P2(:, 2), '-b', P3(:, 1), P3(:, 2), '-r');
p1 = 'see figure 1';

