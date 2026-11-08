clear all; 
close all;
clc;
format long;

name = 'Nathan Phaymany';
id = '20016227';
hw_num = 5;


%% Plotting parameters:
dis_scale_factor = 10;           % amplification factor for displacement 
force_scale_factor = 0.5;           % amplification factor for force
font_size = 16;                    % font size 

%% Set up Engineering Data
E = 210E9;                        % Young's modulus [Pa]
A = 0.25E-4;                            % cross-sectional area [m^2]
DOF = 2;                            % degree of freedom per node
P = 10E3;                            % applied load [N]

%% Set up Geometry, Boundary Conditions and Loading:
% Nodes:
n = 1;
node(n).pos = [0, 0];                 % x-, y- coordinates
node(n).dis =  [0,  0];                 % x-, y- comp. of displacement: NaN = unknown
node(n).force = [NaN, NaN];                % x-, y- comp. of force: NaN = unknown

% Botton Nodes
for n = 2:6
    node(n).pos = [0.05*(n-1), 0];
    node(n).dis = [NaN, NaN];
    node(n).force = [0, -P];
end

% Top Nodes
for n = 12:-1:8
    node(n).pos = [0.05*(13-n), 0.05*3^0.5];
    node(n).dis = [NaN, NaN];
    node(n).force = [0, 0];
end

n = 7;
node(n).pos =  [0.3, 0];             
node(n).dis =  [0, 0];       
node(n).force = [NaN, NaN];

% Elements

% Exterior Elements
for n = 1:11
    element(n).node = [n, n+1];
    element(n).E = E;
    element(n).A = A;
end

% Interior Elements
element_nodes = dictionary( ...
    12, {[1, 12]}, ...
    13, {[2, 12]}, ...
    14, {[3, 12]}, ...
    15, {[3, 11]}, ...
    16, {[3, 10]}, ...
    17, {[4, 10]}, ...
    18, {[5, 10]}, ...
    19, {[5, 9]}, ...
    20, {[5, 8]}, ...
    21, {[6, 8]} );
for n = 1:numel(keys(element_nodes))
    index = 11+n;
    element(index).node = element_nodes{index};
    element(index).E = E;
    element(index).A = A;
end

node_original = node;
element_original = element;

%% Set up FEA model

function [newNode, newElement] = FEA(node, element, DOF)
    % Define a function to compute stiffness matrix for truss element
    K_e_fun = @(E,A,L,beta) E*A/L .* ...
        [cos(beta)^2, cos(beta)*sin(beta), -cos(beta)^2, -cos(beta)*sin(beta);
        cos(beta)*sin(beta), sin(beta)^2, -cos(beta)*sin(beta), -sin(beta)^2;
        -cos(beta)^2, -cos(beta)*sin(beta), cos(beta)^2, cos(beta)*sin(beta);
        -cos(beta)*sin(beta), -sin(beta)^2, -cos(beta)*sin(beta), sin(beta)^2];

    for n = 1:numel(element)
        dx = node(element(n).node(2)).pos(1) - node(element(n).node(1)).pos(1);
        dy = node(element(n).node(2)).pos(2) - node(element(n).node(1)).pos(2);
        % Compute length L
        element(n).L = sqrt(dx^2 + dy^2);
        
        % Obtain slope angle beta with respect to positive x-axis
        element(n).beta = atan2(dy, dx);

        % Evaluate elemental stiffness matrix
        element(n).K_e = K_e_fun(element(n).E, element(n).A, element(n).L, element(n).beta);
                    
        % Expand element stiffness to global dimension
        K_global_e = zeros(DOF * numel(node));
        index = [2*element(n).node(1) - 1, 2*element(n).node(2) - 1];
        K_global_e(index(1):index(1)+1, index(1):index(1)+1) = element(n).K_e(1:2,1:2);
        K_global_e(index(2):index(2)+1, index(2):index(2)+1) = element(n).K_e(3:4,3:4);
        K_global_e(index(1):index(1)+1, index(2):index(2)+1) = element(n).K_e(1:2,3:4);
        K_global_e(index(2):index(2)+1, index(1):index(1)+1) = element(n).K_e(3:4,1:2);

        element(n).K_global_e = K_global_e;
    end

    % Assemble global stiffness matrix
    K_global = zeros( DOF * numel(node) );
    for n = 1:numel(element)
        K_global = K_global + element(n).K_global_e;
    end

    % Apply boundary conditions to reduce the equilibrium equations
    % Assemble global displacement and force vector
    U_global = zeros(DOF * numel(node), 1);
    F_global = zeros(DOF * numel(node), 1);
    for n = 1:numel(node)
        U_global(2*n-1:2*n) = transpose(node(n).dis);
        F_global(2*n-1:2*n) = transpose(node(n).force);
    end

    % Identify row and column to keep: 
    % where the displacement is unknown (i.e. set to NaN during the setup)
    row_col_to_keep = find(isnan(U_global));
    K_reduced = K_global(row_col_to_keep, row_col_to_keep);
    F_reduced = F_global(row_col_to_keep, 1);

    %% Solve equilibrium equations

    U_reduced = K_reduced \ F_reduced;
    U_global(row_col_to_keep) = U_reduced;
    F_global = K_global * U_global;

    %% Postprocessing 
    % Pass updated global displacement and force back to nodes
    for n = 1:numel(node)
        node(n).dis = U_global(2*n-1:2*n)';
        node(n).force = F_global(2*n-1:2*n)';
    end
        
    % Compute axial stress for each element
    for n = 1:numel(element)
        U_element = zeros(4, 1);
        U_element(1:2) = U_global(2*element(n).node(1)-1:2*element(n).node(1));
        U_element(3:4) = U_global(2*element(n).node(2)-1:2*element(n).node(2));
        element(n).stress = element(n).E/element(n).L .* ...
        [-cos(element(n).beta) -sin(element(n).beta) cos(element(n).beta) sin(element(n).beta)] * U_element;
    end

    newNode = node;
    newElement = element;
end

[node, element] = FEA(node_original, element_original, DOF);

%% Figure 1: Plot orignal and deformed trusses 
figure(1);
hold on;

% Element
for n = 1:numel(element)
    x_elem = [node(element(n).node(1)).pos(1), node(element(n).node(2)).pos(1)];
    y_elem = [node(element(n).node(1)).pos(2), node(element(n).node(2)).pos(2)];
    plot(x_elem,y_elem,'--ko','MarkerFaceColor','k','LineWidth',3);

    dx_elem = [node(element(n).node(1)).dis(1), node(element(n).node(2)).dis(1)];
    dy_elem = [node(element(n).node(1)).dis(2), node(element(n).node(2)).dis(2)];
    p1 = patch(x_elem + dis_scale_factor * dx_elem, ...
               y_elem + dis_scale_factor * dy_elem,...
               sqrt(dx_elem.^2 + dy_elem.^2),...
               'EdgeColor','interp','LineWidth',3);
    text(mean(x_elem),mean(y_elem), sprintf('[%d]',n),...
               'Color','b','FontSize',font_size);
end
colorbar; colormap('jet');

% Node and force vectors
for n = 1:numel(node)
    xd = node(n).pos(1);
    yd = node(n).pos(2);
    max_force = max(abs([node(:).force]));
    q = quiver(xd,yd,...
        node(n).force(1)/max_force,node(n).force(2)/max_force, force_scale_factor,...
        'Color','r','LineWidth',2,'MaxHeadSize', force_scale_factor);
    text(xd, yd, sprintf('%d',n),'Color','b','FontSize',font_size);
end
box on; grid on; axis equal;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Color shows displacement in meters (disp. scaling factor of %.0f)',dis_scale_factor));
set(gca,'FontSize',font_size);

p1a = node(4).dis
p1b = element(1).stress
p2c = 'See figure 1'

% Part D
i = 1;
for A = 0.05E-4:0.05E-4:0.5E-4
    node = node_original;
    element = element_original;

    for n = 1:numel(element)
        element(n).A = A;
    end

    [node, element] = FEA(node, element, DOF);
    maxStress(i) = max(abs([element(:).stress]));
    i = i + 1;
end

figure(2);
plot(0.05:0.05:0.5, maxStress*1E-9, '-o');
xlabel('Area A [cm^2]');
ylabel('Maximum Stress [GPa]');
title('Maximum Stress vs Area');
grid on;

p1d = 'See figure 2'