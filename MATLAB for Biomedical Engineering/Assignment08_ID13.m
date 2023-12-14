%% Section 1
preamble
% Define symbolic variable and ODE
syms t y(t) k M
ode1 = k * y * (1 - y/M) == diff(y, t);
% Find and display the general solution
yGeneral = dsolve(ode1);
disp('General Solution(s):')
disp(yGeneral(1))
% Find, display, and plot the particular solution
M1 = 12500;
k1 = 0.04;
yParticular = dsolve(ode1, 'y(0) = 600', 'IgnoreAnalyticConstraints', true);
yParticular = subs(yParticular, [k , M], [k1, M1]);
disp('Particular Solution:')
disp(yParticular)
% Plotting the particular solution
figure;
fplot(yParticular, [0, 300], 'LineWidth', 2);
hold on;
% Set q and compute t1
q = 95;
t1 = solve(subs(yParticular, t, 't1') == q/100 * M1, 'Real', true);
disp(['Time t1 at which the population reaches ', num2str(q), '% of carrying capacity:'])
disp(double(t1))
% Plot vertical line at t1
plot([t1, t1], [0, M1], '--r', 'LineWidth', 1);
% Compute t2 and check if it's a maximum
d2ydt2 = diff(diff(yParticular, t), t);
t2Candidates = solve(d2ydt2 == 0, 'Real', true);
% Select the candidate where the second derivative is negative (local maximum)
validCandidates = t2Candidates(double(subs(d2ydt2, t, t2Candidates)) < 0);
t2 = validCandidates(1);
disp('Time t2 at which the growth rate attains its maximum:')
disp(double(t2))
% Plot vertical line at t2
plot([t2, t2], [0, M1], '--b', 'LineWidth', 1);
% Add axis labels, title, and legend
xlabel('Time (t)');
ylabel('Population Size');
title('Verhulst-Pearl Model');
legend('y1(t)', 't1 - 95\% of carrying capacity', 't2 - Max Growth Rate');
hold off;
%%
%% Section 2
preamble
syms r t x(t)
% Option setting
opt = odeset('MaxStep', 0.002);
% Parameters
b = 3;
s = 10;
r_it = [10; 20; 40];
% Compute r-independent fixed points and critical temperatur
fixed_points_x1 = [0; 0; 0];
critical_temp = s.*(s+b+3)./(s-b-1)
% Iteration of r
for i=1:3
    r_current = r_it(i);
    % Compute r-dependent fixed points
    fixed_points_x2 = [sqrt(b.*(r_current-1)); sqrt(b.*(r_current-1)); r_current-1];
    fixed_points_x3 = [-sqrt(b.*(r_current-1)); -sqrt(b.*(r_current-1)); r_current-1];
    % Compute solution with initial conditions
    [time,result] = ode23(@(t, x) ode(x, s, r_current, b), [0 50], [0.1 0.1 0.1], opt);
    % Plot1
    figure;
    hold on;
    plot3(result(:, 1),result(:, 2),result(:, 3))
    plot3(fixed_points_x1(1), fixed_points_x1(2), fixed_points_x1(3), '*', 'LineWidth', 1)
    plot3(fixed_points_x2(1), fixed_points_x2(2), fixed_points_x2(3), '*', 'LineWidth', 1)
    plot3(fixed_points_x3(1), fixed_points_x3(2), fixed_points_x3(3), '*', 'LineWidth', 1)
    plot3(result(1, 1), result(1, 2), result(1, 3), '+', 'LineWidth', 1)
    plot3(result(end, 1), result(end, 2), result(end, 3), 'o', 'LineWidth', 1)
    set(gca, 'cameraPosition', [mean(xlim), mean(ylim), max(zlim)])
    xlabel('$x$')
    ylabel('$y$')
    zlabel('$z$')
    grid on
    legend({'trajectory', 'fixed point $x_1$', 'fixed point $x_2$', 'fixed point $x_3$', 'first point of trajectory', 'last point of trajectory'})
    title(sprintf('Trajectory and fixed points for $x_0=y_0=z_0=0.1$ and $r=%d$', r_current))
    hold off;
    % Plot2
    figure;
    hold on;
    subplot(3, 1, 1);
    plot(time, result(:, 1))
    xlabel('$t$')
    ylabel('$x$')
    grid on
    title(sprintf('The 3 coordinates plotted seperately versus $t$ for $x_0=y_0=z_0=0.1$ and $r=%d$', r_current))
    legend({'$x$ depenent on $t$'})
    subplot(3, 1, 2);
    plot(time, result(:, 2))
    xlabel('$t$')
    ylabel('$y$')
    grid on
    legend({'$y$ depenent on $t$'})
    subplot(3, 1, 3);
    plot(time, result(:, 3))
    xlabel('$t$')
    ylabel('$z$')
    grid on
    legend({'$z$ depenent on $t$'})
    hold off;
end
% Function for ODE
function F = ode(x, s, r, b)
    F = [-s.*x(1)+s.* x(2);
        r.* x(1)-x(2)-x(1).*x(3);
        -b.*x(3)+x(1).*x(2)];
end














