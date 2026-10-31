clear all; 
close all;
clc;
format long;

name = 'Nathan Phaymany';
id = '20016227';
hw_num = 4;

%% Problem 1

%% Problem 2

L = 0.01;
E = 69E9;
P = 10E6;
A0 = 0.02;
A = @(x) A0 * (1 - x/(2*L));

% Part A 
A1 = 3*A0/4;
l = L;
k = E*A1/L * [1 -1; -1 1];
F2 = P;
u2= F2 / k(2, 2);
p2a = u2

% Part B
A1 = 7*A0/8;
l1 = L/2;
A2 = 5*A0/8;
l2 = L/2;
K1 = E*A1/l1 * [1 -1 0; -1 1 0; 0 0 0];
K2 = E*A2/l2 * [0 0 0; 0 1 -1; 0 -1 1];
K = K1 + K2;
% Excludes first row and column due to BC: u1 = 0
K = K(2:3, 2:3);
F = [0; P];
u = rref([K F]); 
u = u(:, 3);
p2b = u(end)

% Part C
figure(2);
hold on;
grid on;

utrue = 0;
relErr = zeros(20, 1);
for n = 20:-1:1
    K = zeros(n+1, n+1);
    for i = 1:n
        Ai = A((L*(i*2-1))/(2*n)); 
        l = L/n;
        Ki = E*Ai/l * [1 -1; -1 1];
        K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + Ki;
    end

    F = [zeros(n-1, 1); P];
    K = K(2:n+1, 2:n+1);
    u = rref([K F]);
    u = u(:, n+1);
    displacement = u(end);

    if n == 20
        utrue = displacement;
        continue
    end

    relErr(n) = abs((displacement - utrue) / utrue) * 100;
end

plot(1:20, relErr, '-ob');
title('Relative Error for Number of Elements');
xlabel('Elements (n)');
ylabel('Relative Error (%)');
p2c = 'See figure 2';

% Part D
n = find(relErr < 0.01, 1);
p2d = n