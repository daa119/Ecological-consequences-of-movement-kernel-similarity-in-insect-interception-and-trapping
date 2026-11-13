clc; clear; close all;
rng(1);
% -------------------------------
% Parameters
N = 10^2;
R = 40;                  % radius of circular domain (m)
max_steps = 8640;
sub_steps = 100;
trap_radius = 0.05;      % pitfall trap radius (m)

% -------------------------------
% Ecological movement parameters
mean_speed = 2.5;        % mean movement speed (m/h)
dt = 10/60;              % time step (h); 10 min = 1/6 h

% Mean step length over one time step
mean_step = mean_speed * dt;

% For a 2D Gaussian random walk:
% E[l] = sigma * sqrt(pi/2)
sigma = mean_step / sqrt(pi/2);

% Optimal matching probability
p = 0.658;
zeta = sqrt(-2*log(p))*cot((pi/2)*(1-p));
gamma = zeta*sigma;

fprintf('Mean movement speed = %.3f m/h\n', mean_speed);
fprintf('Mean step length = %.3f m\n', mean_step);
fprintf('Gaussian scale sigma = %.3f m\n', sigma);
fprintf('Cauchy scale gamma (p*=%.3f) = %.3f m\n', p, gamma);

% -------------------------------
% Snapshot times
steps_per_day = 144;

snap_steps = [0, ...
              30*steps_per_day, ...
              60*steps_per_day];

snap_titles = { ...
    '{\itt} = 0 days', ...
    '{\itt} = 30 days', ...
    '{\itt} = 60 days'};

%% ==========================================================
% REGULAR 3 x 3 GRID
% ===========================================================

d = 20;

trap_x = -d:d:d;
trap_y = -d:d:d;

[trap_X, trap_Y] = meshgrid(trap_x, trap_y);

trap_centers_grid = ...
    [trap_X(:), trap_Y(:)];

%% ==========================================================
% LINEAR TRANSECT OF 9 TRAPS
% ===========================================================

num_traps = 9;

transect_x = linspace(-32,32,num_traps);
transect_y = zeros(1,num_traps);

trap_centers_transect = ...
    [transect_x(:), transect_y(:)];

%% ==========================================================
% INITIAL POSITIONS
% ===========================================================

theta = 2*pi*rand(N,1);
radii = R*sqrt(rand(N,1));

x = radii.*cos(theta);
y = radii.*sin(theta);

positions_all = [x y];

%% ==========================================================
% REMOVE INDIVIDUALS INITIALLY INSIDE GRID TRAPS
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

fprintf('Grid: initially placed = %d, removed inside traps = %d\n', ...
        N, sum(in_trap_grid));

%% ==========================================================
% REMOVE INDIVIDUALS INITIALLY INSIDE TRANSECT TRAPS
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

fprintf('Transect: initially placed = %d, removed inside traps = %d\n', ...
        N, sum(in_trap_transect));

%% ==========================================================
% GRID: SHORT-DISTANCE DISPERSAL (SRW)
% ===========================================================

positions_srw_grid = positions_init_grid;

snapshots_srw_grid = cell(1,3);
snapshots_srw_grid{1} = positions_srw_grid;

for step = 1:max_steps

    step_vectors = ...
        sigma*randn(size(positions_srw_grid));

    new_positions = ...
        positions_srw_grid + step_vectors;

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
            positions_srw_grid(outside,:) + ...
            step_vectors(outside,:);

        dist_from_center(outside) = ...
            sqrt(new_positions(outside,1).^2 + ...
                 new_positions(outside,2).^2);

        outside = dist_from_center > R;

    end

    % --------------------------------
    % Check trap intersections
    still_alive = ...
        true(size(positions_srw_grid,1),1);

    for ss = 1:sub_steps

        frac = ss/sub_steps;

        intermediate_positions = ...
            positions_srw_grid + ...
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

    positions_srw_grid = ...
        new_positions(still_alive,:);

    if step == snap_steps(2)
        snapshots_srw_grid{2} = positions_srw_grid;

    elseif step == snap_steps(3)
        snapshots_srw_grid{3} = positions_srw_grid;
    end

end

%% ==========================================================
% GRID: LONG-DISTANCE DISPERSAL (LW)
% ===========================================================

positions_lw_grid = positions_init_grid;

snapshots_lw_grid = cell(1,3);
snapshots_lw_grid{1} = positions_lw_grid;

for step = 1:max_steps

    step_lengths = ...
        abs(gamma * ...
        tan(pi*(rand(size(positions_lw_grid,1),1)-0.5)));

    angles = ...
        2*pi*rand(size(positions_lw_grid,1),1);

    dx = step_lengths.*cos(angles);
    dy = step_lengths.*sin(angles);

    step_vectors = [dx dy];

    new_positions = ...
        positions_lw_grid + step_vectors;

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
            positions_lw_grid(outside,:) + ...
            step_vectors(outside,:);

        dist_from_center(outside) = ...
            sqrt(new_positions(outside,1).^2 + ...
                 new_positions(outside,2).^2);

        outside = dist_from_center > R;

    end

    % --------------------------------
    % Check trap intersections
    still_alive = ...
        true(size(positions_lw_grid,1),1);

    for ss = 1:sub_steps

        frac = ss/sub_steps;

        intermediate_positions = ...
            positions_lw_grid + ...
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

    positions_lw_grid = ...
        new_positions(still_alive,:);

    if step == snap_steps(2)
        snapshots_lw_grid{2} = positions_lw_grid;

    elseif step == snap_steps(3)
        snapshots_lw_grid{3} = positions_lw_grid;
    end

end

%% ==========================================================
% TRANSECT: SHORT-DISTANCE DISPERSAL (SRW)
% ===========================================================

positions_srw_transect = ...
    positions_init_transect;

snapshots_srw_transect = cell(1,3);
snapshots_srw_transect{1} = ...
    positions_srw_transect;

for step = 1:max_steps

    step_vectors = ...
        sigma*randn(size(positions_srw_transect));

    new_positions = ...
        positions_srw_transect + step_vectors;

    % --------------------------------
    % Confined circular boundary
    dist_from_center = ...
        sqrt(new_positions(:,1).^2 + ...
             new_positions(:,2).^2);

    outside = dist_from_center > R;

    while any(outside)

        step_vectors(outside,:) = ...
            sigma*randn(sum(outside),2);

        new_positions(outside,:) = ...
            positions_srw_transect(outside,:) + ...
            step_vectors(outside,:);

        dist_from_center(outside) = ...
            sqrt(new_positions(outside,1).^2 + ...
                 new_positions(outside,2).^2);

        outside = dist_from_center > R;

    end

    % --------------------------------
    % Check trap intersections
    still_alive = ...
        true(size(positions_srw_transect,1),1);

    for ss = 1:sub_steps

        frac = ss/sub_steps;

        intermediate_positions = ...
            positions_srw_transect + ...
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

    positions_srw_transect = ...
        new_positions(still_alive,:);

    if step == snap_steps(2)
        snapshots_srw_transect{2} = ...
            positions_srw_transect;

    elseif step == snap_steps(3)
        snapshots_srw_transect{3} = ...
            positions_srw_transect;
    end

end

%% ==========================================================
% TRANSECT: LONG-DISTANCE DISPERSAL (LW)
% ===========================================================

positions_lw_transect = ...
    positions_init_transect;

snapshots_lw_transect = cell(1,3);
snapshots_lw_transect{1} = ...
    positions_lw_transect;

for step = 1:max_steps

    step_lengths = ...
        abs(gamma * ...
        tan(pi*(rand(size(positions_lw_transect,1),1)-0.5)));

    angles = ...
        2*pi*rand(size(positions_lw_transect,1),1);

    dx = step_lengths.*cos(angles);
    dy = step_lengths.*sin(angles);

    step_vectors = [dx dy];

    new_positions = ...
        positions_lw_transect + step_vectors;

    % --------------------------------
    % Confined circular boundary
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
            positions_lw_transect(outside,:) + ...
            step_vectors(outside,:);

        dist_from_center(outside) = ...
            sqrt(new_positions(outside,1).^2 + ...
                 new_positions(outside,2).^2);

        outside = dist_from_center > R;

    end

    % --------------------------------
    % Check trap intersections
    still_alive = ...
        true(size(positions_lw_transect,1),1);

    for ss = 1:sub_steps

        frac = ss/sub_steps;

        intermediate_positions = ...
            positions_lw_transect + ...
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

    positions_lw_transect = ...
        new_positions(still_alive,:);

    if step == snap_steps(2)
        snapshots_lw_transect{2} = ...
            positions_lw_transect;

    elseif step == snap_steps(3)
        snapshots_lw_transect{3} = ...
            positions_lw_transect;
    end

end

%% ==========================================================
% PLOT: 2 x 3 POPULATION DISTRIBUTIONS
% ===========================================================

figure('Color','w', ...
       'Position',[100 100 1450 760]);

%% ----------------------------------------------------------
% TOP ROW: REGULAR GRID
% ----------------------------------------------------------

for k = 1:3

    subplot(2,3,k);
    hold on;

    % Circular confined boundary
    theta_boundary = linspace(0,2*pi,500);

    plot(R*cos(theta_boundary), ...
         R*sin(theta_boundary), ...
         '--k', ...
         'LineWidth',1.0);

    % Initial positions: light green
    if k == 1

        scatter(snapshots_srw_grid{k}(:,1), ...
                snapshots_srw_grid{k}(:,2), ...
                18, ...
                [0.4 0.8 0.4], ...
                'filled');

    else

        % SRW walkers: black
        if ~isempty(snapshots_srw_grid{k})

            scatter(snapshots_srw_grid{k}(:,1), ...
                    snapshots_srw_grid{k}(:,2), ...
                    18, ...
                    'k', ...
                    'filled');

        end

        % LW walkers: blue
        if ~isempty(snapshots_lw_grid{k})

            scatter(snapshots_lw_grid{k}(:,1), ...
                    snapshots_lw_grid{k}(:,2), ...
                    18, ...
                    'b', ...
                    'filled');

        end

    end

    % Pitfall traps: hollow circles
    for t = 1:size(trap_centers_grid,1)

        plot(trap_centers_grid(t,1), ...
             trap_centers_grid(t,2), ...
             'ro', ...
             'MarkerSize',6, ...
             'MarkerFaceColor','none', ...
             'LineWidth',1.2);

    end

    % Legend outside top-right panel
    if k == 3

        hInitial = scatter(NaN, NaN, ...
            18, [0.4 0.8 0.4], 'filled');

        hShort = scatter(NaN, NaN, ...
            18, 'k', 'filled');

        hLong = scatter(NaN, NaN, ...
            18, 'b', 'filled');

        % Store original subplot position
        axPos = get(gca,'Position');

        lgd = legend([hInitial hShort hLong], { ...
            'initial position', ...
            'short-distance dispersal', ...
            'long-distance dispersal'}, ...
            'FontName','Arial', ...
            'FontSize',14, ...
            'Box','off');

        % Restore subplot to its original size
        set(gca,'Position',axPos);

        % Position legend manually outside the top-right panel
        set(lgd,'Units','normalized');
        lgdPos = get(lgd,'Position');

        set(lgd,'Position', ...
            [0.79 0.76 lgdPos(3) lgdPos(4)]);

    end

    axis equal;

    xlim([-R R]);
    ylim([-R R]);

    title(snap_titles{k}, ...
        'Interpreter','tex', ...
        'FontName','Arial', ...
        'FontSize',16, ...
        'FontWeight','normal');

    if k == 1

        ylabel('space, {\ity} (m)', ...
            'Interpreter','tex', ...
            'FontName','Arial', ...
            'FontSize',16);

    end

    box on;

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',14, ...
        'LineWidth',1.0, ...
        'TickDir','out');

    hold off;

end

%% ----------------------------------------------------------
% BOTTOM ROW: TRANSECT
% ----------------------------------------------------------

for k = 1:3

    subplot(2,3,k+3);
    hold on;

    % Circular confined boundary
    theta_boundary = linspace(0,2*pi,500);

    plot(R*cos(theta_boundary), ...
         R*sin(theta_boundary), ...
         '--k', ...
         'LineWidth',1.0);

    % Initial positions: light green
    if k == 1

        scatter(snapshots_srw_transect{k}(:,1), ...
                snapshots_srw_transect{k}(:,2), ...
                18, ...
                [0.4 0.8 0.4], ...
                'filled');

    else

        % SRW walkers: black
        if ~isempty(snapshots_srw_transect{k})

            scatter(snapshots_srw_transect{k}(:,1), ...
                    snapshots_srw_transect{k}(:,2), ...
                    18, ...
                    'k', ...
                    'filled');

        end

        % LW walkers: blue
        if ~isempty(snapshots_lw_transect{k})

            scatter(snapshots_lw_transect{k}(:,1), ...
                    snapshots_lw_transect{k}(:,2), ...
                    18, ...
                    'b', ...
                    'filled');

        end

    end

    % Pitfall traps: hollow circles
    for t = 1:size(trap_centers_transect,1)

        plot(trap_centers_transect(t,1), ...
             trap_centers_transect(t,2), ...
             'ro', ...
             'MarkerSize',6, ...
             'MarkerFaceColor','none', ...
             'LineWidth',1.2);

    end

    axis equal;

    xlim([-R R]);
    ylim([-R R]);

    xlabel('space, {\itx} (m)', ...
        'Interpreter','tex', ...
        'FontName','Arial', ...
        'FontSize',16);

    if k == 1

        ylabel('space, {\ity} (m)', ...
            'Interpreter','tex', ...
            'FontName','Arial', ...
            'FontSize',16);

    end

    box on;

    set(gca, ...
        'FontName','Arial', ...
        'FontSize',14, ...
        'LineWidth',1.0, ...
        'TickDir','out');

    hold off;

end