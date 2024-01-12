%% Section 1
preamble
% Define the equation
syms y(x)
Dy = diff(y);
Dy2 = diff(y, 2);
Dy3 = diff(y, 3);
eq = 3*x^4*Dy3 + x^3*Dy2 - 4*x^2*Dy - 5*sqrt(x)*y;
eqlatex = '3x^4y''''''+x^3y''''-4x^2y''-5\sqrt{x}y=0';
% Convert to 1st order system
[eqs, vars] = reduceDifferentialOrder(eq, y(x));
[M, F] = massMatrixForm(eqs, vars);
f = M\F;
ode = odeFunction(f, vars);
x1 = 1; % lower boundary
x2 = 10; % upper boundary
% Initial guess for y, y', and y''
yini = x.^3; % Third-order polynomial as the initial guess
yinit=matlabFunction([-yini;-simplify(diff(yini));-simplify(diff(yini, 2))]);
% Boundary conditions
bv = [-2; 2; -1]; % y(x1), y(x2), y'(x2)
% Residuals at boundaries
bcres = @(y1, y2, bv) [y1(1) - bv(1); y2(1) - bv(2); y2(2) - bv(3)];
% Solve the boundary value problem
xinit = linspace(x1, x2, 100); % Grid for initial guess
solinit = bvpinit(xinit, yinit);
solnum = bvp4c(ode, @(y1, y2) bcres(y1, y2, bv), solinit);
% Evaluate the solution on a dense grid for plotting and comparison
xs = linspace(x1, x2, 500);
[soly, solyp] = deval(solnum, xs);
ys = soly(1, :); % Solution y
ys1 = solyp(1, :); % First derivative of y
ys2 = solyp(2, :); % Second derivative of y
ys3 = solyp(3, :); % Third derivative of y
% Plotting the solution and its first three derivatives
figure;
hold on;
plot(xs, ys);
plot(xs, solyp);
% Substitute solution into the left-hand side of the differential equation
lhs = 3*xs.^4.*ys3 + xs.^3.*ys2 - 4*xs.^2.*ys1 - 5*sqrt(xs).*ys;
plot(xs, lhs, 'k--');
legend('$y$', '$y''$', '$y''''$', '$y''''''$', 'LHS')
xlabel('$x$')
ylabel('$y(x)$, $y''(x)$, $y''''(x)$, $y''''''$, LHS')
title(sprintf('Solution of $%s$', eqlatex), sprintf('$y(%d)=%d$, $y(%d)=%d$, $y''(%d)=%d$', x1, bv(1), x2, bv(2), x2, bv(2)))
grid on
hold off
% Display the infinity norm
disp("Infinity norm: " + norm(lhs, inf));
%%
%% Section 2
preamble
% Define the parameters
L = [2, 5, 10, 20];
q = 0.5;
x1 = 0; % lower boundary
x2 = 5; % upper boundary
xinit = linspace(x1, x2, 100); % Grid for initial guess
xs = linspace(x1, x2, 500);
% Define the ODE
ode = @(x,y,lam)  [y(2); -(lam - 2 * q * cos(2 * x)) * y(1)];
% Initialize array to store eigenvalues
eigenvalues = zeros(4, 1);
% Boundary conditions
bv = [0; 0; 1]; % y(x1), y(x2), y'(x1)
% Residuals at boundaries
bcres = @(y1, y2, lam, bv) [y1(1) - bv(1); y2(1) - bv(2); y1(2) - bv(3)];
% Initial guess
yinit=@(x, lam) [cos(x * sqrt(lam)); -sqrt(lam) * sin(x * sqrt(lam))];
for i = 1:4
    % Solve the boundary value problem
    solinit=bvpinit(xinit, @(x) yinit(x, L(i)), L(i));
    sol = bvp4c(ode, @(y1, y2, lam) bcres(y1, y2, lam, bv), solinit);
    % Extract the eigenvalue
    eigenvalues(i) = double(sol.parameters);
    % Evaluate the eigenfunction and its derivatives
    [soly,solyp] = deval(sol, xs);
    ys = soly(1, :); % solution y
    ys2 = solyp(2, :); % second derivative of y
    % Plot the eigenfunction and its derivatives
    figure;
    hold on;
    plot(xs, ys); 
    plot(xs, solyp);
    % Calculate the left-hand side of the equation and add to the plot
    lhs = ys2 + (eigenvalues(i) - 2 * q * cos(2 .* xs)) .* ys;
    plot(xs, lhs);
    % Check the infinity norm
    assert(norm(lhs, inf) <= 1e-3, 'Infinity norm exceeds 1e-3');
    % Add labels, title, and legend
    xlabel('$x$');
    ylabel('Function Value');
    title(sprintf('Mathieu''s ODE with $\\lambda=%s$', num2str(eigenvalues(i))));
    legend('Eigenfunction', '1st Derivative', '2nd Derivative', 'LHS')
    grid on;
    hold off;
end
% Display the eigenvalues
disp('Eigenvalues:');
disp(eigenvalues);


