%% T1D Glucose-Insulin Simulation
clear;
clc;
close all;

%% Parameters (Subject 1)

b0 = 0.95;
b1 = 0.001;
b2 = 0.0002;

c1 = 0.05;
c2 = 0.04;

rb = 0.30;

%% Equilibrium Point

% Target euglycemic glucose level
G_star = 80;

% Insulin concentration required to maintain G* = 80 mg/dL
I_star = (b0 - b1*G_star) / (b2*G_star);

% Required basal insulin infusion rate in mU/min
rb_mU_min = c1*I_star/c2;

% Convert to U/hr for interpretation
rb_U_hr = rb_mU_min*60/1000;

fprintf('Equilibrium glucose G* = %.4f mg/dL\n', G_star);
fprintf('Equilibrium insulin I* = %.4f mU/L\n', I_star);
fprintf('Required basal infusion = %.4f mU/min\n', rb_mU_min);
fprintf('Required basal infusion = %.4f U/hr\n', rb_U_hr);

%% Jacobian Matrix

J_star = [
    -(b1 + b2*I_star),   -b2*G_star;
    0,                  -c1
    ];

eigenvalues = eig(J_star);

disp('Jacobian matrix:')
disp(J_star)

disp('Eigenvalues:')
disp(eigenvalues)

%% Baseline Simulation

% Initial conditions near the equilibrium point
G0 = 120;
I0 = 65;

% Convert to perturbation variables
g0 = G0 - G_star;
i0 = I0 - I_star;

y0 = [g0; i0];

% Simulation time
tspan = [0 300];

% Linearized model
linear_model = @(t,y) J_star*y;

% Solve the system
[t,y] = ode45(linear_model, tspan, y0);

% Convert back to original variables
G = G_star + y(:,1);
I = I_star + y(:,2);

%% Figure 1: Glucose and Insulin Dynamics

figure;

subplot(2,1,1)

plot(t,G,'LineWidth',2);
hold on;
yline(G_star,'--','Equilibrium G^*','LineWidth',1.5);

xlabel('Time (min)');
ylabel('Glucose (mg/dL)');
title('Plasma Glucose Dynamics Near Equilibrium');
grid on;

subplot(2,1,2)

plot(t,I,'LineWidth',2);
hold on;
yline(I_star,'--','Equilibrium I^*','LineWidth',1.5);

xlabel('Time (min)');
ylabel('Insulin (mU/L)');
title('Plasma Insulin Dynamics Near Equilibrium');
grid on;

%% Figure 2: Different Initial Conditions

% Three initial conditions near the equilibrium point
initial_conditions = [
    90,  58;
    120, 65;
    70,  50
    ];

figure;
hold on;

for k = 1:size(initial_conditions,1)

    G0_current = initial_conditions(k,1);
    I0_current = initial_conditions(k,2);

    % Convert to perturbation variables
    y0_current = [
        G0_current - G_star;
        I0_current - I_star
        ];

    % Solve the linearized system
    [t_current,y_current] = ode45( ...
        linear_model, ...
        tspan, ...
        y0_current ...
        );

    % Convert back to glucose
    G_current = G_star + y_current(:,1);

    % Plot each trajectory
    plot(t_current,G_current, ...
        'LineWidth',2, ...
        'DisplayName', ...
        sprintf('G(0)=%.0f, I(0)=%.0f', ...
        G0_current,I0_current));
end

% Add equilibrium line
yline(G_star,'--','Equilibrium G^*', ...
    'LineWidth',1.5);

xlabel('Time (min)');
ylabel('Glucose (mg/dL)');
title('Glucose Trajectories Under Different Initial Conditions');
legend('Location','best');
grid on;
hold off;

%% Figure 3: Sensitivity Analysis for b2

b2_values = [
    0.5*b2;
    1.0*b2;
    1.5*b2
    ];

figure;
hold on;

for k = 1:length(b2_values)

    b2_current = b2_values(k);

    % Keep the same equilibrium glucose level
    G_star_current = 80;

    % Recalculate equilibrium insulin
    I_star_current = ...
        (b0 - b1*G_star_current) ...
        /(b2_current*G_star_current);

    % Recalculate the Jacobian matrix
    J_current = [
        -(b1 + b2_current*I_star_current), ...
        -b2_current*G_star_current;
        0, ...
        -c1
        ];

    % Use the same relative perturbation
    y0_current = [40; 10];

    model_current = @(t,y) J_current*y;

    [t_current,y_current] = ode45( ...
        model_current, ...
        tspan, ...
        y0_current ...
        );

    G_current = G_star_current + y_current(:,1);

    plot(t_current,G_current, ...
        'LineWidth',2, ...
        'DisplayName', ...
        sprintf('b_2 = %.5f',b2_current));
end

yline(80,'--','Equilibrium G^*', ...
    'LineWidth',1.5);

xlabel('Time (min)');
ylabel('Glucose (mg/dL)');
title('Sensitivity of Glucose Dynamics to b_2');
legend('Location','best');
grid on;
hold off;


%% Figure 4: Sensitivity Analysis for c1

c1_values = [
    0.5*c1;
    1.0*c1;
    1.5*c1
    ];

figure;
hold on;

for k = 1:length(c1_values)

    c1_current = c1_values(k);

    % Keep the same equilibrium glucose level
    G_star_current = 80;

    % Keep the same equilibrium insulin level
    I_star_current = I_star;

    % Recalculate the basal infusion rate required
    rb_current = c1_current*I_star_current/c2;

    % Recalculate the Jacobian matrix
    J_current = [
        -(b1 + b2*I_star_current), ...
        -b2*G_star_current;
        0, ...
        -c1_current
        ];

    % Use the same perturbation from equilibrium
    y0_current = [40; 10];

    model_current = @(t,y) J_current*y;

    [t_current,y_current] = ode45( ...
        model_current, ...
        tspan, ...
        y0_current ...
        );

    % Convert back to original variables
    G_current = G_star_current + y_current(:,1);

    plot(t_current,G_current, ...
        'LineWidth',2, ...
        'DisplayName', ...
        sprintf('c_1 = %.3f',c1_current));
end

yline(G_star,'--','Equilibrium G^*', ...
    'LineWidth',1.5);

xlabel('Time (min)');
ylabel('Glucose (mg/dL)');
title('Sensitivity of Glucose Dynamics to c_1');
legend('Location','best');
grid on;
hold off;

exportgraphics(gcf, ...
    'Figure4_c1Sensitivity.png', ...
    'Resolution',300);
