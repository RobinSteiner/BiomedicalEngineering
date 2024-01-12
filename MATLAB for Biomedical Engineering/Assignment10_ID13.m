preamble
% Load the data
load pardata.mat;
% Objective functions
obj_funcs = {
    @(p) mean((X(:,2) - (p(1)*X(:,1).^2 + p(2)*X(:,1) + p(3))).^2),
    @(p) mean(abs(X(:,2) - (p(1)*X(:,1).^2 + p(2)*X(:,1) + p(3)))),
    @(p) median((X(:,2) - (p(1)*X(:,1).^2 + p(2)*X(:,1) + p(3))).^2)
};
% Initial guess
p0 = [1, 1, 1];
% Minimizers (output is turned off for fminunc)
options = optimoptions('fminunc', 'Display', 'off');
minimizers = {{'fminsearch', @fminsearch}, {'fminunc', @(func, p) fminunc(func, p, options)}};
% Loop over minimizers
for i = 1:length(minimizers)
    fprintf('***** Minimizer is %s\n', minimizers{i}{1})
    p = p0;  % Starting point
    % Plotting
    figure;
    hold on;
    plot(X(:,1), X(:,2), '.');
    % Loop over objective functions
    for j = 1:length(obj_funcs)
        % Minimization
        [p, fval] = minimizers{i}{2}(obj_funcs{j}, p);
        % Display results
        fprintf('(A) Objfun at minimum: %.4f Parabola parameters a,b,c: %.4f  %.4f  %.4f\n', fval, p(1), p(2), p(3));
        % Plotting
        x_fit = linspace(min(X(:,1)), max(X(:,1)), 100);
        y_fit = p(1)*x_fit.^2 + p(2)*x_fit + p(3);
        plot(x_fit, y_fit, '-');
    end
    hold off;
    grid on;
    legend('Data points', 'Parabola (A)', 'Parabola (B)', 'Parabola (C)');
    title(sprintf('Parabolas fitted with %s', minimizers{i}{1}));
    xlabel('$x$'); ylabel('$y$');
end