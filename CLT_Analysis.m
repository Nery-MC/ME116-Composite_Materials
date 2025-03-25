function CLT_Analysis
    % CLT Analysis for Composite Laminate
    
    % Layer thickness
    t = 0.2; % Thickness of each layer (mm)
    num_layers = 8;
    
    % Define cumulative thickness positions (mid-plane centered)
    h = linspace(-num_layers*t/2, num_layers*t/2, num_layers+1);
    
    % Example ply orientations (in degrees)
    angles = [0, 90, 0, 90, 0, 90, 0, 90]; 

    % Material properties (modify as needed)
    E1 = 15.6374; % Young's modulus in GPa (fiber direction)
    E2 = 5.1426;  % Young's modulus in GPa (transverse direction)
    G12 = 1.9728; % Shear modulus in GPa
    v12 = 0.32;    % Poisson's ratio
    
    % Convert angles to radians
    theta = deg2rad(angles);
    
    % Initialize A, B, D matrices
    A = zeros(3, 3);
    B = zeros(3, 3);
    D = zeros(3, 3);
    
    % Initialize strains and curvatures
    strains = zeros(num_layers, 3);
    curvatures = zeros(num_layers, 3);
    
    % Compute stiffness matrices
    for i = 1:num_layers
        Q = calculate_Q_matrix(theta(i), E1, E2, G12, v12);
        A = A + Q * (h(i+1) - h(i));
        B = B + 0.5 * Q * (h(i+1)^2 - h(i)^2);
        D = D + (1/3) * Q * (h(i+1)^3 - h(i)^3);
        
        % Calculate strains and curvatures (assuming unit loading)
        strains(i, :) = Q * [0; 0; 1] * (h(i+1) - h(i));
        curvatures(i, :) = Q * [0; 1; 0] * (h(i+1) - h(i));
    end
    
    % Display results
    disp('A matrix:'); disp(A);
    disp('B matrix:'); disp(B);
    disp('D matrix:'); disp(D);

    % Display strains and curvatures
    disp('Strains:'); disp(strains);
    disp('Curvatures:'); disp(curvatures);

    % Plot stress distribution
    plot_stress_distribution(theta, strains, h, E1, E2, G12, v12);

    % Calculate laminate engineering constants
    calculate_laminate_constants(A, B, D, h);
end

function Q = calculate_Q_matrix(theta, E1, E2, G12, v12)
    % Calculate the reduced stiffness matrix (Q) in the local coordinates

    % Calculate compliance matrix
    v21 = (v12 * E2) / E1;  % Reciprocal Poisson's ratio
    Q11 = E1 / (1 - v12 * v21);
    Q12 = (v12 * E2) / (1 - v12 * v21);
    Q22 = E2 / (1 - v12 * v21);
    Q66 = G12;

    % Reduced stiffness matrix
    Q = [Q11, Q12, 0;
         Q12, Q22, 0;
         0,    0,  Q66];

    % Rotation matrix
    c = cos(theta);
    s = sin(theta);
    T = [c^2, s^2, 2*c*s;
         s^2, c^2, -2*c*s;
         -c*s, c*s, c^2 - s^2];

    % Transform Q to the global coordinate system
    Q = T \ Q / T;
end

function plot_stress_distribution(theta, strains, h, E1, E2, G12, v12)
    % Plot stress distribution through thickness
    num_layers = length(theta);
    stress = zeros(num_layers, 3);

    for i = 1:num_layers
        Q = calculate_Q_matrix(theta(i), E1, E2, G12, v12);
        stress(i, :) = Q * strains(i, :)';
    end

    % Plot
    figure;
    for i = 1:3
        subplot(3, 1, i);
        plot(h(1:end-1), stress(:, i), '-o');
        xlabel('Thickness (mm)');
        ylabel(['Stress ' num2str(i)]);
        title(['Stress Distribution: Component ' num2str(i)]);
    end

    % Display average stress components
    avg_stress_x = mean(stress(:, 1));
    avg_stress_y = mean(stress(:, 2));
    avg_stress_xy = mean(stress(:, 3));

    disp(['Average Stress X: ', num2str(avg_stress_x)]);
    disp(['Average Stress Y: ', num2str(avg_stress_y)]);
    disp(['Average Shear Stress XY: ', num2str(avg_stress_xy)]);
end

function calculate_laminate_constants(A, B, D, h)
    % Calculate laminate engineering constants (Ex, Ey, vxy, Gxy)

    h_total = max(h) - min(h);  % Total laminate thickness

    % Compute in-plane stiffness matrix
    Q_bar = [A, B; B, D] / h_total;

    % Compute laminate engineering constants
    Ex = Q_bar(1, 1);
    Ey = Q_bar(2, 2);
    vxy = Q_bar(1, 2) / Q_bar(1, 1);
    Gxy = Q_bar(3, 3);

    % Display results
    disp('Laminate Engineering Constants:');
    disp(['Ex: ' num2str(Ex) ' GPa']);
    disp(['Ey: ' num2str(Ey) ' GPa']);
    disp(['νxy: ' num2str(vxy)]);
    disp(['Gxy: ' num2str(Gxy) ' GPa']);
end
