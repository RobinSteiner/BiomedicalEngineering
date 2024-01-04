% Read the data from the file
data = dlmread('histogram_Y73.txt', ';', 1, 0);

% Extract the X and Y values from the data
X = data(:, 1);
Y = data(:, 2);

% Histogram plot
figure;
bar(X, Y);
xlabel('GV-Domain +/- 0.5');
ylabel('Counts');
title('Histogram Plot');
grid on;

% Calculate the probability density function (PDF)
PDF = normalize(Y, 'norm');

% Probability density plot
figure;
bar(X, PDF);
xlabel('GV-Domain +/- 0.5');
ylabel('Probability Density');
title('Probability Density Function');
grid on;

[gv_flat, gv_soft, gv_thr, gv_max] = findGVs(PDF);

% Probability density plot with special points
figure;
bar(X, PDF);
xlabel('GV-Domain +/- 0.5');
ylabel('Probability Density');
title('Probability Density Function');
hold on;
plot(gv_flat(1) - 0.5, gv_flat(2), 'bo', 'MarkerSize', 8, 'Marker', 'x', 'Color', 'r');
plot(gv_soft(1) - 0.5, gv_soft(2), 'bo', 'MarkerSize', 8, 'Marker', '+', 'Color', 'r');
plot(gv_thr(1) - 0.5, gv_thr(2), 'bo', 'MarkerSize', 8, 'Marker', '*', 'Color', 'r');
plot(gv_max(1) - 0.5, gv_max(2), 'bo', 'MarkerSize', 8, 'Marker', 'o', 'Color', 'r');
legend('PDF', "GV_{flat} = (" + gv_flat(1) + ", " + gv_flat(2) + ")", ...
"GV_{soft} = (" + gv_soft(1) + ", " + gv_soft(2) + ")", ...
"GV_{thr} = (" + gv_thr(1) + ", " + gv_thr(2) + ")", ...
"GV_{max} = (" + gv_max(1) + ", " + gv_max(2) + ")");
grid on;

% Values from literature
gv_max_minus_one = 156;
l_voxel = 0.324;
l_cort = 0.23;

% Extracting the part between GV_thr and GV_max which is associated with the bone
X_bone = X(gv_thr(1) : gv_max(1));
Y_bone = normalize(Y(gv_thr(1) : gv_max(1)), 'norm');

gv_ev = (l_voxel / l_cort) * (gv_max(1) + gv_max_minus_one - 2 * gv_soft(1)) + gv_soft(1)
phi_vas = (X_bone - gv_ev) / (gv_soft(1) - gv_ev);

% Plotting vascular porosity
figure;
plot(100 * phi_vas, Y_bone);
hold on
grid on
title('Probability density function of vascular porosity')
ylabel('Probability')
xlabel('Porosity in [%]')

% Calculating longitudinal young's modulus
hom_matrices = arrayfun(@(phi) inv(hom_exvas_to_macro(phi)), phi_vas, 'UniformOutput', false);
stiff_matrix_33 = cellfun(@(mat) 1/mat(3,3), hom_matrices);

% Plotting longitudinal young's modulus
figure;
bar(stiff_matrix_33, Y_bone)
grid on
title('Youngs modulus')
xlabel('Longitudinal Young´s modulus [GPa]')
ylabel('number of voxels')
grid on


function[gv_flat, gv_soft, gv_thr, gv_max] = findGVs(PDF)
    threshold = 0.001;
    maxima = 1;
    minima = 1;
    for i = 50:(length(PDF) - 1)
        if PDF(i) > PDF(i - 1) && PDF(i) > PDF(i + 1) && PDF(i) - PDF(i - 3) > threshold
            if(maxima == 1)
                gv_flat = [i; PDF(i)];
            elseif (maxima == 2)
                gv_soft = [i; PDF(i)];
            end
            maxima = maxima + 1;
        elseif PDF(i) < PDF(i - 1) && PDF(i) < PDF(i + 1) && PDF(i - 3) - PDF(i) > threshold
            if(minima == 2)
                gv_thr = [i; PDF(i)];
            end
            minima = minima + 1;
        elseif PDF(i) == 0
            gv_max = [i; PDF(i)];
            break;
        end
    end

    disp("GV_flat =")
    disp(gv_flat)
    disp("GV_soft =")
    disp(gv_soft)
    disp("GV_thr =")
    disp(gv_thr)
    disp("GV_max =")
    disp(gv_max)
end


