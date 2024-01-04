% Define the parameters
L = [2, 5, 10, 20];
q = 0.5;
x1 = 1; % lower boundary
x2 = 10; % upper boundary
xinit = linspace(0, 5, 100);
xs = linspace(0, 5, 500);
% Define the ODE
mat = @(x,y,lam)  [y(2); -(lam - 2 * q * cos(2 * x)) * y(1)];
% Initialize array to store eigenvalues
eigenvalues = zeros(4, 1);
% Boundary conditions
bv = [0; 0; 1]; % y(x1), y(x2), y'(x1)
% Residuals at boundaries
bcres = @(y1, y2, lam, bv) [y1(1) - bv(1); y2(1) - bv(2); y2(2) - bv(3)];
% Initial guess
yinit=@(x, lam) [cos(x * sqrt(lam)); -sqrt(lam)*sin(x*sqrt(lam))];
for i = 1:4
    % Solve the boundary value problem
    solinit=bvpinit(xinit, @(x) yinit(x, L(i)), L(i));
    sol = bvp4c(mat, @(y1, y2, lam) bcres(y1, y2, lam, bv), solinit);
    % Extract the eigenvalue
    eigenvalues(i) = double(sol.parameters);
    % Evaluate the eigenfunction and its derivatives
    [soly,solyp] = deval(sol, xs);
    ys=soly(1,:); % solution y
    ys2=solyp(2,:); % second derivative of y
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

