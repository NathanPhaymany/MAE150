clear all; 
close all;
clc;
format long;

name = 'Nathan Phaymany';
id = '20016227';
hw_num = 2;

%% Plotting Function

histPlot = @(data) histogram(data, 100, 'Normalization', 'pdf');

%% Problem 1

rng('default');
n = 25000;

figure;
% Lift Force (uniform)
ax1 = subplot(2, 6, [1, 2]);
pd1 = makedist('Uniform', 'lower', 134.5, 'upper', 135.5);
L = random(pd1, [1,n]);
histPlot(L);

% Fluid density (triangular)
ax2 = subplot(2, 6, [3, 4]);
pd2 = makedist('Triangular', 'a', 1.189, 'b', 1.19, 'c', 1.191);
rho = random(pd2, [1,n]);
histPlot(rho);

% Airspeed (triangular)
ax3 = subplot(2, 6, [5, 6]);
pd3 = makedist('Triangular', 'a', 32.9, 'b', 33.4, 'c', 33.9);
v = random(pd3, [1,n]);
histPlot(v);

% Air Foil Area (gaussian)
ax4 = subplot(2, 6, [7, 9]);
pd4 = makedist('normal', 'mu', 5, 'sigma', 0.05);
S = random(pd4, [1,n]);
histPlot(S);

% Lift Coefficient
ax5 = subplot(2, 6, [10, 12]);
C = L ./ (0.5 * rho .* v.^2 .* S);
histPlot(C);

ax = [ax1, ax2, ax3, ax4, ax5];
label = ["Lift Force (L)", "Fluid Density (\rho)", "Airspeed (v)", "Air Foil Area (S)", "Lift Coefficient (C_L)"];
for i = 1:length(ax)
    xlabel(ax(i), label(i));
    ylabel(ax(i), 'PDF');
    title(ax(i), label(i));
    box(ax(i), 'on');
    grid(ax(i), 'on');
end

p1a = 'See figure 1'

phat = mle(C);
p1b = phat(1)
p1c = phat(2)
p1d = normcdf(0.05, p1b, p1c)


%% Problem 2

% files = ['screw.STL', 'handle.STL', 'brace.STL', "wheel.STL"];

% for i = 1:length(files)
%     model = stlread(files(i));
%     figure;
%     trimesh(model, 'FaceColor', 'none', 'EdgeColor', 'k');
%     axis equal;
% end

% p2 = 'See figure 2'