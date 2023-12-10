%% Section 1
preamble
% Define symbolic functions
syms x y
f = 1 - x^2;
g = 0.5 - (x^2)/3;
% Convert symbolic expressions to functions
f_num = matlabFunction(f);
g_num = matlabFunction(g);
% Plot functions
figure;
fplot(f, [-1, 1], 'DisplayName', 'f(x)');
hold on;
fplot(g, [-1, 1], 'DisplayName', 'g(x)');
legend('Location', 'southwest');
xlabel('x');
ylabel('y');
title('Intersection of Functions');
% Compute intersection points
intersection_points = solve(f == g, x);
% Fill the region R
fill_region_x = linspace(double(intersection_points(1)), double(intersection_points(2)), 100);
fill_region_y1 = 1 - fill_region_x.^2;
fill_region_y2 = 0.5 - (fill_region_x.^2)/3;
fill([fill_region_x, fliplr(fill_region_x)], [fill_region_y1, fliplr(fill_region_y2)], 'c', 'FaceAlpha', 0.2, 'DisplayName', 'Region R');
% Compute and display area by symbolic integration
area_symbolic = int(g - f, x, intersection_points(1), intersection_points(2));
disp('Symbolic Area:');
disp(char(area_symbolic) + " ~= " + double(area_symbolic));
% Compute and display area by numerical integration
area_numerical = integral(@(x) g_num(x) - f_num(x), double(intersection_points(1)), double(intersection_points(2)));
disp('Numerical Area:');
disp(num2str(area_numerical));
% Calculate relative error
relative_error_area = abs((area_numerical - double(area_symbolic)) / double(area_symbolic));
disp('Relative Error of Area:');
disp([num2str(relative_error_area), newline]);
% Compute and display centroid by symbolic integration
centroid_x_symbolic = int(x*(g - f), x, intersection_points(1), intersection_points(2))/area_symbolic;
centroid_y_symbolic = int((g + f)*(g - f)/2, x, intersection_points(1), intersection_points(2))/area_symbolic;
disp('Symbolic Centroid:');
disp("x: " + char(centroid_x_symbolic) + " ~= " + double(centroid_x_symbolic));
disp("y: " + char(centroid_y_symbolic) + " ~= " + double(centroid_y_symbolic));
% Compute and display centroid by numerical integration
centroid_x_numerical = integral(@(x) x.*(g_num(x) - f_num(x)), double(intersection_points(1)), double(intersection_points(2)))/area_numerical;
centroid_y_numerical = integral(@(x) (g_num(x) + f_num(x)).*(g_num(x) - f_num(x))/2, double(intersection_points(1)), double(intersection_points(2)))/area_numerical;
disp('Numerical Centroid:');
disp(['x: ', num2str(centroid_x_numerical)]);
disp(['y: ', num2str(centroid_y_numerical), newline]);
% Plot the centroid
plot(centroid_x_numerical, centroid_y_numerical, 'ro', 'MarkerSize', 10, 'DisplayName', 'Numerical Centroid');
plot(double(centroid_x_symbolic), double(centroid_y_symbolic), 'bx', 'MarkerSize', 10, 'DisplayName', 'Symbolic Centroid');
legend('Location', 'southwest');
% Define function for h(x, y)
h = x^2 + y^2;
% Compute and display integral of h over region R by numerical integration
integral_h_numeric = integral2(matlabFunction(h), double(intersection_points(1)), double(intersection_points(2)), f_num, g_num);
disp('Numerical Integral of h:');
disp(num2str(integral_h_numeric));
% Compute and display integral of h over region R by symbolic integration
integral_h_symbolic = int(int(h, y, f, g), x, intersection_points(1), intersection_points(2));
disp('Symbolic Integral of h:');
disp(char(integral_h_symbolic) + " ~= " + double(integral_h_symbolic));
% Calculate relative error
relative_error_integral_h = abs((integral_h_numeric - double(integral_h_symbolic)) / double(integral_h_symbolic));
disp('Relative Error of Numerical h Integral:');
disp([num2str(relative_error_integral_h), newline]);
% Define function for arc length
arc_length_f = sqrt(1 + diff(f, x)^2);
arc_length_g = sqrt(1 + diff(g, x)^2);
% Compute and display circumference by symbolic integration
circumference_symbolic = int(arc_length_f, x, intersection_points(1), intersection_points(2)) + ...
    int(arc_length_g, x, intersection_points(1), intersection_points(2));
disp('Symbolic Circumference:');
disp(char(circumference_symbolic) + " ~= " + double(circumference_symbolic));
% Compute and display circumference by numerical integration
circumference_numerical = integral(matlabFunction(arc_length_f), double(intersection_points(1)), double(intersection_points(2))) + ...
    integral(matlabFunction(arc_length_g), double(intersection_points(1)), double(intersection_points(2)));
disp('Numerical Circumference:');
disp(num2str(circumference_numerical));
% Calculate relative error
relative_error_circumference = abs((circumference_numerical - double(circumference_symbolic)) / double(circumference_symbolic));
disp('Relative Error of Numerical Circumference:');
disp(num2str(relative_error_circumference));
%% Section 2
preamble
% Define symbolic variables
syms a t real positive
% Define symbolic functions
x = 2*a*(1 - cos(t))*cos(t);
y = 2*a*(1 - cos(t))*sin(t);
% Plot the cardioid with a = 2
a_value = 2;
figure;
fplot(subs(x, a, a_value), subs(y, a, a_value), [0, 2*pi], 'DisplayName', 'Cardioid with a = 2');
title('Cardioid for a = 2');
xlabel('x');
ylabel('y');
axis equal;
legend;
% Compute dx/dt and dy/dt
dx_dt = diff(x, t);
dy_dt = diff(y, t);
% Compute the differential of the path length ds
ds = sqrt(dx_dt^2 + dy_dt^2);
% Compute and display curve length of C by symbolic integration
curve_length_symbolic = int(ds, t, 0, 2*pi);
disp('Symbolic Curve Length:');
disp(curve_length_symbolic);
% Compute and display curve length of C by numerical integration with a=2
curve_length_numerical = integral(matlabFunction(subs(ds, a, 2)), 0, 2*pi);
disp('Numerical Curve Length:');
disp(curve_length_numerical);
% Compute and display area enclosed by C by symbolic integration
area = y * dx_dt;
area_symbolic = int(area, t, 2*pi, 0);
disp('Symbolic Area:');
disp(area_symbolic);
% Compute and display area enclosed by C by numerical integration with a=2
area_numerical = integral(matlabFunction(subs(area, a, 2)), 2*pi, 0);
disp('Numerical Area:');
disp(area_numerical);
