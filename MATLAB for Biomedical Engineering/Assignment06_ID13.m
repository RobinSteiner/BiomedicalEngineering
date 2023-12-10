%% Section 1
preamble
% Define symbolic variables
syms x y real
% Define the equations for the curve C and the ellipse E
C = (x^2 - x*y/2 + y^2)^2 - 2*(x^2 - y^2);
E = 4*x^2 + y^2 - 1;
% Find the symbolic intersection points
intersectionPoints = solve([C, E], [x, y]);
% Extract the coordinates of the intersection point in the second quadrant
intersectionPoints = double([intersectionPoints.x intersectionPoints.y]);
intersectionPoint = intersectionPoints(intersectionPoints(:, 1) < 0 & intersectionPoints(:, 2) > 0, :);
% Calculate the slopes of the tangents at the intersection point
slopeC = subs(-diff(C, x) / diff(C, y), [x y], intersectionPoint);
slopeE = subs(-diff(E, x) / diff(E, y), [x y], intersectionPoint);
% Calculate the angles between the tangents
alpha = atan(abs((slopeC - slopeE) / (1 + slopeC * slopeE)));
angleDegrees = double(rad2deg(alpha));
% Plot the curves and tangents and the angle
figure;
hold on;
fimplicit(C, [-2, 2, -2, 2], 'b', 'LineWidth', 2);
fimplicit(E, [-2, 2, -2, 2], 'r', 'LineWidth', 2);
refline(double(slopeC), -double(slopeC) * intersectionPoint(1) + intersectionPoint(2));
refline(double(slopeE), -double(slopeE) * intersectionPoint(1) + intersectionPoint(2));
scatter(intersectionPoint(1), intersectionPoint(2), 50, 'k', 'filled');
text(intersectionPoint(1) - 0.5, intersectionPoint(2) + 0.3, "$\alpha = " + num2str(angleDegrees, 4) + "^\circ$", 'FontSize', 12);
hold off;
% Set axis labels, title, and legend
xlabel('x');
ylabel('y');
title('Intersection of C and E with Tangents');
legend({'C', 'E', 'Tangent to C', 'Tangent to E', 'Intersection Point'}, 'Location', 'NorthWest');
% Enforce the same scale and range on both axes
axis equal;
axis([-2, 2, -2, 2]);