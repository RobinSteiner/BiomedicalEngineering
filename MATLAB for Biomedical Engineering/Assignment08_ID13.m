preamble
% Define symbolic variable and ODE
syms t y(t) k M
ode = k * y * (1 - y/M) == diff(y, t);
% Find and display the general solution
yGeneral = dsolve(ode);
disp('General Solution(s):')
disp(yGeneral(1))
% Find, display, and plot the particular solution
M1 = 12500;
k1 = 0.04;
yParticular = dsolve(ode, 'y(0) = 600', 'IgnoreAnalyticConstraints', true);
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
