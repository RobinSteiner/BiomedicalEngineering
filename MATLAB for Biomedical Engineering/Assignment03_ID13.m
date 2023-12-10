%% Section 1
preamble
X = 30;
delta = 1e-4;
f = @(x) cosint(x) .* cos(x);
figure;
x_values = linspace(0, X, 1000);
plot(x_values, f(x_values), 'b');
hold on;
% Find zeros of f(x) in the interval [0, X]
starting_points = 0.5:0.5:X;
zeros = arrayfun(@(xs) fzero(@(x) f(x), xs), starting_points);
rounded_zeros = round(zeros, 12);
unique_zeros = unique(rounded_zeros);
plot(unique_zeros, f(unique_zeros), 'ro', 'MarkerSize', 8);
% Find local minima and maxima
xmin = unique_zeros(1:end-1);
xmax = unique_zeros(2:end);
local_minima = arrayfun(@(a, b) fminbnd(@(x) f(x), a, b), xmin, xmax);
local_minima = local_minima(abs(f(local_minima)) > delta);
plot(local_minima, f(local_minima), 'g+', 'MarkerSize', 8);
local_maxima = arrayfun(@(a, b) fminbnd(@(x) -f(x), a, b), xmin, xmax);
local_maxima = local_maxima(abs(f(local_maxima)) > delta);
plot(local_maxima, f(local_maxima), 'b+', 'MarkerSize', 8);
% Add grid, x-axis, title, axis labels, and legend
grid on;
xlabel('$x$');
ylabel('$f(x)$');
ylim([-0.3, 0.3])
title('Plot of $f(x) = cosint(x) \cdot cos(x)$');
legend('$f(x)$', 'Zeros', 'Local Minima', 'Local Maxima');
hold off;