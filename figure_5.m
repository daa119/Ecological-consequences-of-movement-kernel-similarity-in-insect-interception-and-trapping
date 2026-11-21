% Mean daily cumulative pitfall capture over 30 days for a
% grid and a transect. Each of 10 realisations uses the same starting draw
% for both layouts before trap-specific exclusions. Each capture percentage
% uses its own post-exclusion population. Shading shows pointwise 95%
% Student t intervals across the realisations, bounded by 0% and 100%.
%
% Gaussian displacements have independent components and Rayleigh
% step lengths. The long-distance case uses folded-Cauchy step lengths
% and uniform headings. Proposals ending outside the circular domain are
% redrawn. An individual is captured if any of the 100 sampled locations
% on the accepted segment lies in a trap, including a check at the endpoint.
% Short intersections between sampled locations may go undetected.

clc; clear; close all;
rng(25);  % Fixed seed for this figure.

save_figure = true;

% -------------------------------
% Parameters
N = 10^2;
num_runs = 10;           % use at least two runs for the confidence bands
R = 40;                  % radius of circular domain (m)
max_steps = 4320;
sub_steps = 100;         % equally spaced checks per accepted segment
trap_radius = 0.05;      % pitfall trap radius (m)

% -------------------------------
% Ecological movement parameters
mean_speed = 2.5;        % mean movement speed (m/h)
dt = 10/60;              % time step (h); 10 min = 1/6 h

% Mean step length over one time step
mean_step = mean_speed * dt;

% For a 2D SRW with independent Gaussian increments,
% the radial step length is Rayleigh distributed with
% E[l] = sigma * sqrt(pi/2)
sigma = mean_step / sqrt(pi/2);

fprintf('Mean movement speed = %.3f m/h\n', mean_speed);
fprintf('Time step = %.3f h\n', dt);
fprintf('Mean step length = %.3f m\n', mean_step);
fprintf('Gaussian scale sigma = %.3f m\n', sigma);

% -------------------------------
% Match tail probabilities. The p = 0.658 case minimises the L^2 distance.

p1 = 0.1;
zeta1 = sqrt(-2*log(p1))*cot((pi/2)*(1-p1));
gamma1 = zeta1*sigma;

p2 = 0.5;
zeta2 = sqrt(-2*log(p2))*cot((pi/2)*(1-p2));
gamma2 = zeta2*sigma;

p3 = 0.658;              % optimal survival probability
zeta3 = sqrt(-2*log(p3))*cot((pi/2)*(1-p3));
gamma3 = zeta3*sigma;

p4 = 0.9;
zeta4 = sqrt(-2*log(p4))*cot((pi/2)*(1-p4));
gamma4 = zeta4*sigma;

gamma_values = [gamma1, gamma2, gamma3, gamma4];

fprintf('gamma (p = 0.1)   = %.3f m\n', gamma1);
fprintf('gamma (p = 0.5)   = %.3f m\n', gamma2);
fprintf('gamma (p*= 0.658) = %.3f m\n', gamma3);
fprintf('gamma (p = 0.9)   = %.3f m\n', gamma4);

%% ==========================================================
% PANEL (a): REGULAR 3 x 3 GRID
% ===========================================================

d = 20;                  % inter-trap spacing (m)

trap_x = -d:d:d;
trap_y = -d:d:d;

[trap_X, trap_Y] = meshgrid(trap_x, trap_y);

trap_centers_grid = ...
    [trap_X(:), trap_Y(:)];

%% ==========================================================
% PANEL (b): LINEAR TRANSECT OF 9 TRAPS
% ===========================================================

num_traps = 9;

transect_x = linspace(-32,32,num_traps);
transect_y = zeros(1,num_traps);

trap_centers_transect = ...
    [transect_x(:), transect_y(:)];

%% ==========================================================
% STORAGE ACROSS SIMULATION RUNS
% ===========================================================

steps_per_day = 144;
num_days = 30;
days = 1:num_days;

trap_counts_gauss_grid_runs = zeros(num_days,num_runs);
trap_counts_levy_grid_runs = zeros(num_days,length(gamma_values),num_runs);
trap_counts_gauss_transect_runs = zeros(num_days,num_runs);
trap_counts_levy_transect_runs = zeros(num_days,length(gamma_values),num_runs);

initial_population_grid_runs = zeros(num_runs,1);
initial_population_transect_runs = zeros(num_runs,1);
initial_removed_grid_runs = zeros(num_runs,1);
initial_removed_transect_runs = zeros(num_runs,1);

for run = 1:num_runs

%% ==========================================================
% INITIAL POSITIONS
% ===========================================================

% Square-root radius gives uniform density per unit area.
theta = 2*pi*rand(N,1);
radii = R*sqrt(rand(N,1));

x = radii.*cos(theta);
y = radii.*sin(theta);

positions_all = [x y];

%% ==========================================================
% PANEL (a): REMOVE INDIVIDUALS INITIALLY INSIDE GRID TRAPS
% ===========================================================

in_trap_grid = false(N,1);

for i = 1:size(trap_centers_grid,1)

    dist_to_trap = ...
        sqrt((positions_all(:,1)-trap_centers_grid(i,1)).^2 + ...
             (positions_all(:,2)-trap_centers_grid(i,2)).^2);

    in_trap_grid = ...
        in_trap_grid | (dist_to_trap < trap_radius);

end

positions_init_grid = ...
    positions_all(~in_trap_grid,:);

initial_population_grid = ...
    size(positions_init_grid,1);

initial_population_grid_runs(run) = initial_population_grid;
initial_removed_grid_runs(run) = sum(in_trap_grid);

fprintf('Run %d, Grid: initially placed = %d, removed inside traps = %d, effective N0 = %d\n', ...
        run, N, sum(in_trap_grid), initial_population_grid);

%% ==========================================================
% PANEL (a): GAUSSIAN RANDOM WALK
% ===========================================================

positions_gauss_grid = positions_init_grid;
trapped_gauss_grid = zeros(max_steps,1);

for step = 1:max_steps

    step_vectors = ...
        sigma*randn(size(positions_gauss_grid));

    new_positions = ...
        positions_gauss_grid + step_vectors;

    % --------------------------------
    % Confined circular boundary
    % Redraw any step ending outside R
    dist_from_center = ...
        sqrt(new_positions(:,1).^2 + ...
             new_positions(:,2).^2);

    outside = dist_from_center > R;

    while any(outside)

        step_vectors(outside,:) = ...
            sigma*randn(sum(outside),2);

        new_positions(outside,:) = ...
            positions_gauss_grid(outside,:) + ...
            step_vectors(outside,:);

        dist_from_center(outside) = ...
            sqrt(new_positions(outside,1).^2 + ...
                 new_positions(outside,2).^2);

        outside = dist_from_center > R;

    end

    % --------------------------------
    % Check the accepted segment for trap contact at the sampled points.
    still_alive = ...
        true(size(positions_gauss_grid,1),1);

    for ss = 1:sub_steps

        frac = ss/sub_steps;

        intermediate_positions = ...
            positions_gauss_grid + ...
            frac*step_vectors;

        for t = 1:size(trap_centers_grid,1)

            dist_to_trap = ...
                sqrt((intermediate_positions(:,1) - ...
                      trap_centers_grid(t,1)).^2 + ...
                     (intermediate_positions(:,2) - ...
                      trap_centers_grid(t,2)).^2);

            still_alive = ...
                still_alive & ...
                (dist_to_trap > trap_radius);

        end

    end

    positions_gauss_grid = ...
        new_positions(still_alive,:);

    trapped_gauss_grid(step) = ...
        1 - (size(positions_gauss_grid,1) / ...
        initial_population_grid);

end

%% ==========================================================
% PANEL (a): LEVY WALKS
% ===========================================================

trapped_levy_grid = ...
    zeros(max_steps,length(gamma_values));

for idx = 1:length(gamma_values)

    gamma = gamma_values(idx);
    positions_levy = positions_init_grid;

    for step = 1:max_steps

        step_lengths = ...
            abs(gamma * ...
            tan(pi*(rand(size(positions_levy,1),1)-0.5)));

        angles = ...
            2*pi*rand(size(positions_levy,1),1);

        dx = step_lengths.*cos(angles);
        dy = step_lengths.*sin(angles);

        step_vectors = [dx dy];

        new_positions = ...
            positions_levy + step_vectors;

        % --------------------------------
        % Confined circular boundary
        % Redraw any step ending outside R
        dist_from_center = ...
            sqrt(new_positions(:,1).^2 + ...
                 new_positions(:,2).^2);

        outside = dist_from_center > R;

        while any(outside)

            new_lengths = ...
                abs(gamma * ...
                tan(pi*(rand(sum(outside),1)-0.5)));

            new_angles = ...
                2*pi*rand(sum(outside),1);

            step_vectors(outside,1) = ...
                new_lengths.*cos(new_angles);

            step_vectors(outside,2) = ...
                new_lengths.*sin(new_angles);

            new_positions(outside,:) = ...
                positions_levy(outside,:) + ...
                step_vectors(outside,:);

            dist_from_center(outside) = ...
                sqrt(new_positions(outside,1).^2 + ...
                     new_positions(outside,2).^2);

            outside = dist_from_center > R;

        end

        % --------------------------------
        % Check the accepted segment for trap contact at the sampled points.
        still_alive = ...
            true(size(positions_levy,1),1);

        for ss = 1:sub_steps

            frac = ss/sub_steps;

            intermediate_positions = ...
                positions_levy + ...
                frac*step_vectors;

            for t = 1:size(trap_centers_grid,1)

                dist_to_trap = ...
                    sqrt((intermediate_positions(:,1) - ...
                          trap_centers_grid(t,1)).^2 + ...
                         (intermediate_positions(:,2) - ...
                          trap_centers_grid(t,2)).^2);

                still_alive = ...
                    still_alive & ...
                    (dist_to_trap > trap_radius);

            end

        end

        positions_levy = ...
            new_positions(still_alive,:);

        trapped_levy_grid(step,idx) = ...
            1 - (size(positions_levy,1) / ...
            initial_population_grid);

    end

end

%% ==========================================================
% PANEL (b): REMOVE INDIVIDUALS INITIALLY INSIDE TRANSECT TRAPS
% ===========================================================

in_trap_transect = false(N,1);

for i = 1:size(trap_centers_transect,1)

    dist_to_trap = ...
        sqrt((positions_all(:,1)-trap_centers_transect(i,1)).^2 + ...
             (positions_all(:,2)-trap_centers_transect(i,2)).^2);

    in_trap_transect = ...
        in_trap_transect | (dist_to_trap < trap_radius);

end

positions_init_transect = ...
    positions_all(~in_trap_transect,:);

initial_population_transect = ...
    size(positions_init_transect,1);

initial_population_transect_runs(run) = initial_population_transect;
initial_removed_transect_runs(run) = sum(in_trap_transect);

fprintf('Run %d, Transect: initially placed = %d, removed inside traps = %d, effective N0 = %d\n', ...
        run, N, sum(in_trap_transect), initial_population_transect);

%% ==========================================================
% PANEL (b): GAUSSIAN RANDOM WALK
% ===========================================================

positions_gauss_transect = ...
    positions_init_transect;

trapped_gauss_transect = ...
    zeros(max_steps,1);

for step = 1:max_steps

    step_vectors = ...
        sigma*randn(size(positions_gauss_transect));

    new_positions = ...
        positions_gauss_transect + step_vectors;

    % --------------------------------
    % Confined circular boundary
    % Redraw any step ending outside R
    dist_from_center = ...
        sqrt(new_positions(:,1).^2 + ...
             new_positions(:,2).^2);

    outside = dist_from_center > R;

    while any(outside)

        step_vectors(outside,:) = ...
            sigma*randn(sum(outside),2);

        new_positions(outside,:) = ...
            positions_gauss_transect(outside,:) + ...
            step_vectors(outside,:);

        dist_from_center(outside) = ...
            sqrt(new_positions(outside,1).^2 + ...
                 new_positions(outside,2).^2);

        outside = dist_from_center > R;

    end

    % --------------------------------
    % Check the accepted segment for trap contact at the sampled points.
    still_alive = ...
        true(size(positions_gauss_transect,1),1);

    for ss = 1:sub_steps

        frac = ss/sub_steps;

        intermediate_positions = ...
            positions_gauss_transect + ...
            frac*step_vectors;

        for t = 1:size(trap_centers_transect,1)

            dist_to_trap = ...
                sqrt((intermediate_positions(:,1) - ...
                      trap_centers_transect(t,1)).^2 + ...
                     (intermediate_positions(:,2) - ...
                      trap_centers_transect(t,2)).^2);

            still_alive = ...
                still_alive & ...
                (dist_to_trap > trap_radius);

        end

    end

    positions_gauss_transect = ...
        new_positions(still_alive,:);

    trapped_gauss_transect(step) = ...
        1 - (size(positions_gauss_transect,1) / ...
        initial_population_transect);

end

%% ==========================================================
% PANEL (b): LEVY WALKS
% ===========================================================

trapped_levy_transect = ...
    zeros(max_steps,length(gamma_values));

for idx = 1:length(gamma_values)

    gamma = gamma_values(idx);
    positions_levy = positions_init_transect;

    for step = 1:max_steps

        step_lengths = ...
            abs(gamma * ...
            tan(pi*(rand(size(positions_levy,1),1)-0.5)));

        angles = ...
            2*pi*rand(size(positions_levy,1),1);

        dx = step_lengths.*cos(angles);
        dy = step_lengths.*sin(angles);

        step_vectors = [dx dy];

        new_positions = ...
            positions_levy + step_vectors;

        % --------------------------------
        % Confined circular boundary
        % Redraw any step ending outside R
        dist_from_center = ...
            sqrt(new_positions(:,1).^2 + ...
                 new_positions(:,2).^2);

        outside = dist_from_center > R;

        while any(outside)

            new_lengths = ...
                abs(gamma * ...
                tan(pi*(rand(sum(outside),1)-0.5)));

            new_angles = ...
                2*pi*rand(sum(outside),1);

            step_vectors(outside,1) = ...
                new_lengths.*cos(new_angles);

            step_vectors(outside,2) = ...
                new_lengths.*sin(new_angles);

            new_positions(outside,:) = ...
                positions_levy(outside,:) + ...
                step_vectors(outside,:);

            dist_from_center(outside) = ...
                sqrt(new_positions(outside,1).^2 + ...
                     new_positions(outside,2).^2);

            outside = dist_from_center > R;

        end

        % --------------------------------
        % Check the accepted segment for trap contact at the sampled points.
        still_alive = ...
            true(size(positions_levy,1),1);

        for ss = 1:sub_steps

            frac = ss/sub_steps;

            intermediate_positions = ...
                positions_levy + ...
                frac*step_vectors;

            for t = 1:size(trap_centers_transect,1)

                dist_to_trap = ...
                    sqrt((intermediate_positions(:,1) - ...
                          trap_centers_transect(t,1)).^2 + ...
                         (intermediate_positions(:,2) - ...
                          trap_centers_transect(t,2)).^2);

                still_alive = ...
                    still_alive & ...
                    (dist_to_trap > trap_radius);

            end

        end

        positions_levy = ...
            new_positions(still_alive,:);

        trapped_levy_transect(step,idx) = ...
            1 - (size(positions_levy,1) / ...
            initial_population_transect);

    end

end

%% ==========================================================
% CUMULATIVE TRAP CAPTURE AT EACH DAY
% ===========================================================

trap_counts_gauss_grid = zeros(num_days,1);

trap_counts_levy_grid = ...
    zeros(num_days,length(gamma_values));

trap_counts_gauss_transect = zeros(num_days,1);

trap_counts_levy_transect = ...
    zeros(num_days,length(gamma_values));

for day = 1:num_days

    idx_step = day*steps_per_day;

    trap_counts_gauss_grid(day) = ...
        trapped_gauss_grid(idx_step)*100;

    trap_counts_gauss_transect(day) = ...
        trapped_gauss_transect(idx_step)*100;

    for j = 1:length(gamma_values)

        trap_counts_levy_grid(day,j) = ...
            trapped_levy_grid(idx_step,j)*100;

        trap_counts_levy_transect(day,j) = ...
            trapped_levy_transect(idx_step,j)*100;

    end

end


    trap_counts_gauss_grid_runs(:,run) = trap_counts_gauss_grid;
    trap_counts_levy_grid_runs(:,:,run) = trap_counts_levy_grid;
    trap_counts_gauss_transect_runs(:,run) = trap_counts_gauss_transect;
    trap_counts_levy_transect_runs(:,:,run) = trap_counts_levy_transect;

end

%% ==========================================================
% MEAN CUMULATIVE CAPTURE ACROSS SIMULATION RUNS
% ===========================================================

trap_counts_gauss_grid = mean(trap_counts_gauss_grid_runs,2);
trap_counts_levy_grid = mean(trap_counts_levy_grid_runs,3);
trap_counts_gauss_transect = mean(trap_counts_gauss_transect_runs,2);
trap_counts_levy_transect = mean(trap_counts_levy_transect_runs,3);

%% ==========================================================
% PLOT
% ===========================================================


%% ==========================================================
% PRINT CUMULATIVE CAPTURE AT t = 30 DAYS
% ===========================================================

day30 = num_days;

fprintf('\nMean cumulative trap capture at t = 30 days across %d runs:\n', num_runs);

fprintf('Grid, short-distance: %.1f%%\n', ...
        trap_counts_gauss_grid(day30));

fprintf('Grid, p = 0.1:        %.1f%%\n', ...
        trap_counts_levy_grid(day30,1));

fprintf('Grid, p = 0.5:        %.1f%%\n', ...
        trap_counts_levy_grid(day30,2));

fprintf('Grid, p*= 0.658:      %.1f%%\n', ...
        trap_counts_levy_grid(day30,3));

fprintf('Grid, p = 0.9:        %.1f%%\n', ...
        trap_counts_levy_grid(day30,4));

fprintf('Transect, short-distance: %.1f%%\n', ...
        trap_counts_gauss_transect(day30));

fprintf('Transect, p = 0.1:        %.1f%%\n', ...
        trap_counts_levy_transect(day30,1));

fprintf('Transect, p = 0.5:        %.1f%%\n', ...
        trap_counts_levy_transect(day30,2));

fprintf('Transect, p*= 0.658:      %.1f%%\n', ...
        trap_counts_levy_transect(day30,3));

fprintf('Transect, p = 0.9:        %.1f%%\n', ...
        trap_counts_levy_transect(day30,4));

%% ==========================================================
% POINTWISE 95% CONFIDENCE INTERVALS ACROSS SIMULATION RUNS
% ===========================================================

% Two-sided 95% Student t critical value for the number of runs.
% Invert the beta representation of the Student t tail in base MATLAB.
df = num_runs - 1;
beta_quantile = betaincinv(0.05,df/2,0.5);
t95 = sqrt(df * (1/beta_quantile - 1));

se_gauss_grid = std(trap_counts_gauss_grid_runs,0,2)/sqrt(num_runs);
se_levy_grid = std(trap_counts_levy_grid_runs,0,3)/sqrt(num_runs);
se_gauss_transect = std(trap_counts_gauss_transect_runs,0,2)/sqrt(num_runs);
se_levy_transect = std(trap_counts_levy_transect_runs,0,3)/sqrt(num_runs);

CI_low_grid = [ ...
    trap_counts_gauss_grid - t95*se_gauss_grid, ...
    trap_counts_levy_grid - t95*se_levy_grid];

CI_high_grid = [ ...
    trap_counts_gauss_grid + t95*se_gauss_grid, ...
    trap_counts_levy_grid + t95*se_levy_grid];

CI_low_transect = [ ...
    trap_counts_gauss_transect - t95*se_gauss_transect, ...
    trap_counts_levy_transect - t95*se_levy_transect];

CI_high_transect = [ ...
    trap_counts_gauss_transect + t95*se_gauss_transect, ...
    trap_counts_levy_transect + t95*se_levy_transect];

CI_low_grid(CI_low_grid < 0) = 0;
CI_high_grid(CI_high_grid > 100) = 100;
CI_low_transect(CI_low_transect < 0) = 0;
CI_high_transect(CI_high_transect > 100) = 100;

fig = figure('Color','w', ...
       'Position',[100 100 1200 500]);

%% ---------------- PANEL (a): REGULAR GRID ----------------

subplot(1,2,1);
hold on;

% Pointwise 95% confidence regions across simulation runs
fill([days fliplr(days)], ...
     [CI_low_grid(:,1)' fliplr(CI_high_grid(:,1)')], ...
     [0.5 0.5 0.5], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.10, ...
     'HandleVisibility','off');

for j = 2:5
    if j == 4
        band_color = [1 0 0];
    else
        band_color = [0 0 1];
    end

    fill([days fliplr(days)], ...
         [CI_low_grid(:,j)' fliplr(CI_high_grid(:,j)')], ...
         band_color, ...
         'EdgeColor','none', ...
         'FaceAlpha',0.08, ...
         'HandleVisibility','off');
end

h1 = plot(days, trap_counts_gauss_grid, '-', ...
    'Color','k', ...
    'LineWidth',1.8);

h2 = plot(days, trap_counts_levy_grid(:,1), '--', ...
    'Color','b', ...
    'LineWidth',1.8);

h3 = plot(days, trap_counts_levy_grid(:,2), '-.', ...
    'Color','b', ...
    'LineWidth',1.8);

h4 = plot(days, trap_counts_levy_grid(:,3), '-', ...
    'Color','r', ...
    'LineWidth',1.8);

h5 = plot(days, trap_counts_levy_grid(:,4), ':', ...
    'Color','b', ...
    'LineWidth',1.8);

% Dummy invisible handles used as legend headings
hHeaderShort = plot(NaN,NaN, ...
    'LineStyle','none', ...
    'Marker','none', ...
    'Color','w');

hHeaderLong = plot(NaN,NaN, ...
    'LineStyle','none', ...
    'Marker','none', ...
    'Color','w');

xlabel('time (days)', ...
    'FontSize',16, ...
    'FontName','Arial');

ylabel('cumulative trap capture, {\it C}({\itt}) (%)', ...
    'Interpreter','tex', ...
    'FontSize',16, ...
    'FontName','Arial');

title('(a) regular grid', ...
    'FontSize',16, ...
    'FontName','Arial', ...
    'FontWeight','normal');

legend([hHeaderShort h1 hHeaderLong h2 h3 h4 h5], { ...
    'short-distance movement', ...
    sprintf('   \\sigma = %.3f', sigma), ...
    'long-distance movement', ...
    sprintf('   \\color{blue}{\\itp} = %.1f', p1), ...
    sprintf('   \\color{blue}{\\itp} = %.1f', p2), ...
    sprintf('   \\color{red}{\\itp}^{*} = %.3f', p3), ...
    sprintf('   \\color{blue}{\\itp} = %.1f', p4)}, ...
    'Location','northwest', ...
    'Interpreter','tex', ...
    'FontSize',14, ...
    'FontName','Arial', ...
    'Box','off');

box on;

xlim([0 num_days+1]);
ylim([0 100]);

set(gca,'XTick',[1 5 10 15 20 25 30]);

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out');

hold off;

%% ---------------- PANEL (b): TRANSECT ----------------

subplot(1,2,2);
hold on;

% Pointwise 95% confidence regions across simulation runs
fill([days fliplr(days)], ...
     [CI_low_transect(:,1)' fliplr(CI_high_transect(:,1)')], ...
     [0.5 0.5 0.5], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.10, ...
     'HandleVisibility','off');

for j = 2:5
    if j == 4
        band_color = [1 0 0];
    else
        band_color = [0 0 1];
    end

    fill([days fliplr(days)], ...
         [CI_low_transect(:,j)' fliplr(CI_high_transect(:,j)')], ...
         band_color, ...
         'EdgeColor','none', ...
         'FaceAlpha',0.08, ...
         'HandleVisibility','off');
end

plot(days, trap_counts_gauss_transect, '-', ...
    'Color','k', ...
    'LineWidth',1.8);

plot(days, trap_counts_levy_transect(:,1), '--', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(days, trap_counts_levy_transect(:,2), '-.', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(days, trap_counts_levy_transect(:,3), '-', ...
    'Color','r', ...
    'LineWidth',1.8);

plot(days, trap_counts_levy_transect(:,4), ':', ...
    'Color','b', ...
    'LineWidth',1.8);

xlabel('time (days)', ...
    'FontSize',16, ...
    'FontName','Arial');

title('(b) transect', ...
    'FontSize',16, ...
    'FontName','Arial', ...
    'FontWeight','normal');

box on;

xlim([0 num_days+1]);
ylim([0 100]);

set(gca,'XTick',[1 5 10 15 20 25 30]);

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out');

hold off;

%% Export the PNG with a small white margin around all plotted content
if save_figure
    png_file = fullfile(fileparts(mfilename('fullpath')), 'figure_5.png');
    set(fig, 'PaperPositionMode', 'auto');
    dpi = 300;
    drawnow;
    print(fig, png_file, '-dpng', sprintf('-r%d',dpi));

    rgb = imread(png_file);
    content = any(rgb < 255, 3);
    rows = find(any(content, 2));
    cols = find(any(content, 1));
    if ~isempty(rows) && ~isempty(cols)
        pad = 12;
        r = max(1, rows(1)-pad):min(size(rgb, 1), rows(end)+pad);
        c = max(1, cols(1)-pad):min(size(rgb, 2), cols(end)+pad);
        pixels_per_metre = round(dpi / 0.0254);
        imwrite(rgb(r,c,:), png_file, 'png', ...
            'ResolutionUnit','meter','XResolution',pixels_per_metre);
    end
end
