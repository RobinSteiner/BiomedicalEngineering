%% Boundary value problem with Euler-Cauchy ODE
preamble

%% 3x^4y''' + x^3y'' − 4x^2y' − 5sqrt(x)y = 0, y(1) = -2, y(10) = 2, y'(10) = -1
close all

% Define the equation
syms y(x)
Dy = diff(y);
Dy2 = diff(y, 2);
Dy3 = diff(y, 3);

eq = 3*x^4*Dy3 + x^3*Dy2 - 4*x^2*Dy - 5*sqrt(x)*y;
eqlatex = "$3x^4y'''+x^3y''-4x^2y'-5\sqrt{x}y=0$";

% Convert to 1st order system
[eqs, vars] = reduceDifferentialOrder(eq, y(x));
[M, F] = massMatrixForm(eqs, vars);
f = M\F;
fh = odeFunction(f, vars);
x1 = 1; % lower boundary
x2 = 10; % upper boundary

% Initial guess for y, y', and y''
yinit = @(x) [x.^3; -2*x; -6]; % Third-order polynomial as the initial guess

% Boundary conditions
bv = [-2; 2; -1]; % y(x1), y(x2), y'(x2)
bcond = [y(x1) == bv(1), y(x2) == bv(2), Dy(x2) == bv(3)];

% Residuals at boundaries
bcres = @(y1, y2, y3, bv) [y1(1) - bv(1); y2(1) - bv(2); y2(2) - bv(3)];

titstr = ["Solution of ", eqlatex, "$y(", num2str(x1), ")=", num2str(bv(1)), ",\ y(", num2str(x2), ")=", num2str(bv(2)), ",\ y'(", num2str(x2), ")=", num2str(bv(3)), "$"];

% Solve numerically
xinit = linspace(x1, x2, 100); % Grid for initial guess
solinit = bvpinit(xinit, yinit);
solnum = bvp4c(fh, @(y1, y2, y3) bcres(y1, y2, y3, bv), solinit); % Numerical solver

% Evaluate the solution on a dense grid for plotting and comparison
xs = linspace(x1, x2, 500);
[soly, solyp] = deval(solnum, xs);
ys = soly(1, :); % Solution y
ys1 = soly(2, :); % First derivative of y
ys2 = solyp(2, :); % Second derivative of y

% Plot results
figure
hold on
hnum = plot(xs, ys, 'b');
plot(xs, ys1, 'g'); % Plot y'
plot(xs, ys2, 'r'); % Plot y''
xlabel('$x$')
ylabel('$y(x)$, $y''(x)$, $y''''(x)$')
legend('Numerical solution', '$y''$', '$y''''$')
title(join(titstr), 'fontsize', 16)
grid on
box on

% Substitute solution into the left-hand side of the differential equation
lhs = 3*xs.^4.*ys2 + xs.^3.*ys1 - 4*xs.^2.*ys - 5*sqrt(xs).*ys;
hproof = plot(xs, lhs, 'k--');

% Display the infinity norm
inf_norm = max(abs(lhs));
disp("Maximal absolute difference num-sym: " + inf_norm);

keyboard % Transfer control to keyboard
