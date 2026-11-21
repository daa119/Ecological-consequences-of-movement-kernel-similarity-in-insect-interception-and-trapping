% Spatial distributions at days 0, 10 and 20 for Gaussian movement
% and its optimally linked Cauchy case. Both start from the same population.
% Bars give the percentage of the initial population per 2.5 m bin.
% The printed summaries count emigrants beyond the plotted window.
%
% The starting location determines the movement rule. Cash-crop insects
% are biased toward the trap crop. Within the trap crop both cases use
% the same reduced-scale, unbiased Gaussian walk. Outside the field
% insects use their original kernel at full scale, with unbiased headings
% and possible re-entry. A cash-to-trap crossing triggers an entry-removal
% trial even if the step crosses the entire trap strip. Surviving endpoints
% inside the trap crop receive a separate within-trap removal trial.

clear;
clc;

%% ==============================================================
%% Reproducibility
%% ==============================================================

rng(1);  % Reproducible random sequence within each figure.

save_figure = true;


%% ==============================================================
%% Parameters
%% ==============================================================

dt      = 0.5;        % time step (h)
Npart   = 10000;       % initial population
p_steps = 1441;       % 30 days at dt = 0.5 h

sigma = 0.5;          % Gaussian scale (m)
L     = 100;          % total field length (m)

trap_frac = 0.80;     % cash-crop fraction
zone_left = trap_frac * L;

%% ==============================================================
%% Attraction
%% ==============================================================

beta = 0.1;           % maximum attraction strength
nu   = 1;             % linear spatial attraction profile

%% ==============================================================
%% Arrestment in trap crop
%% ==============================================================

rho = 1/3;

sigma_trap = rho * sigma;

%% ==============================================================
%% Removal
%% ==============================================================

alpha_entry = 0.1;

alpha_day = 0.2;

alpha_step = ...
    1 - (1 - alpha_day)^(dt/24);

%% ==============================================================
%% Gaussian-Cauchy relationship
%% ==============================================================

p_opt = 0.721;

gamma = sqrt(2) * erfcinv(p_opt) / ...
        tan((pi/2)*(1-p_opt)) * sigma;

fprintf('Optimal p = %.3f\n',p_opt);
fprintf('sigma = %.3f m\n',sigma);
fprintf('gamma = %.3f m\n',gamma);
fprintf('gamma/sigma = %.3f\n',gamma/sigma);

%% ==============================================================
%% Snapshot times
%% ==============================================================

snap_days = [0 10 20];

snap_steps = ...
    round(snap_days*24/dt) + 1;

%% ==============================================================
%% Common initial condition
%% ==============================================================

x0 = zone_left * rand(Npart,1);

%% ==============================================================
%% SHORT-DISTANCE MOVEMENT: GAUSSIAN
%% ==============================================================

x = x0;

is_killed = false(Npart,1);

x_snap_gauss = cell(1,length(snap_steps));
x_snap_gauss{1} = x;

for a = 2:p_steps

    active = ~is_killed;

    %% Position before movement

    x_old = x;

    was_in_cash = ...
        active & ...
        (x_old >= 0) & ...
        (x_old <= zone_left);

    in_trap = ...
        active & ...
        (x_old > zone_left) & ...
        (x_old <= L);

    outside = ...
        active & ...
        ((x_old < 0) | (x_old > L));

    %% Step magnitude

    step_mag = zeros(Npart,1);

    % Cash crop: Gaussian movement
    if any(was_in_cash)

        step_mag(was_in_cash) = ...
            abs(randn(sum(was_in_cash),1)) .* sigma;

    end

    % Trap crop: common arrested Gaussian movement
    if any(in_trap)

        step_mag(in_trap) = ...
            abs(randn(sum(in_trap),1)) .* sigma_trap;

    end

    % Outside field: Gaussian movement resumes
    if any(outside)

        step_mag(outside) = ...
            abs(randn(sum(outside),1)) .* sigma;

    end

    %% Direction

    dir = zeros(Npart,1);

    % Attraction in cash crop

    x_clipped = ...
        max(0,min(zone_left,x_old));

    beta_space = ...
        (x_clipped ./ zone_left).^nu;

    beta_eff = ...
        beta .* beta_space;

    toward_prob = ...
        (1 + beta_eff)/2;

    r = rand(Npart,1);

    dir(was_in_cash) = -1;

    dir(was_in_cash & (r < toward_prob)) = 1;

    % Unbiased direction in trap crop and outside field

    unbiased_idx = ...
        in_trap | outside;

    if any(unbiased_idx)

        r2 = rand(sum(unbiased_idx),1);

        temp_dir = ones(sum(unbiased_idx),1);

        temp_dir(r2 > 0.5) = -1;

        dir(unbiased_idx) = temp_dir;

    end

    %% Position update

    x_new = x_old;

    x_new(active) = ...
        x_old(active) + ...
        dir(active).*step_mag(active);

    %% Entry-triggered removal

    crossed_into_trap = ...
        was_in_cash & ...
        (x_new > zone_left);

    to_kill_entry = ...
        crossed_into_trap & ...
        (rand(Npart,1) < alpha_entry);

    x = x_new;

    if any(to_kill_entry)

        is_killed(to_kill_entry) = true;
        x(to_kill_entry) = NaN;

    end

    %% Within-trap removal

    active = ~is_killed;

    in_strip = ...
        active & ...
        (x > zone_left) & ...
        (x <= L);

    to_kill_trap = ...
        in_strip & ...
        (rand(Npart,1) < alpha_step);

    if any(to_kill_trap)

        is_killed(to_kill_trap) = true;
        x(to_kill_trap) = NaN;

    end

    %% Save snapshots

    snap_idx = find(snap_steps == a);

    if ~isempty(snap_idx)

        x_snap_gauss{snap_idx} = x;

    end

end

%% ==============================================================
%% LONG-DISTANCE MOVEMENT: CAUCHY
%% ==============================================================

x = x0;

is_killed = false(Npart,1);

x_snap_levy = cell(1,length(snap_steps));
x_snap_levy{1} = x;

for a = 2:p_steps

    active = ~is_killed;

    %% Position before movement

    x_old = x;

    was_in_cash = ...
        active & ...
        (x_old >= 0) & ...
        (x_old <= zone_left);

    in_trap = ...
        active & ...
        (x_old > zone_left) & ...
        (x_old <= L);

    outside = ...
        active & ...
        ((x_old < 0) | (x_old > L));

    %% Step magnitude

    step_mag = zeros(Npart,1);

    % Cash crop: Cauchy movement
    if any(was_in_cash)

        U = rand(sum(was_in_cash),1);

        step_mag(was_in_cash) = ...
            gamma .* tan((pi/2).*U);

    end

    % Trap crop: common arrested Gaussian movement
    if any(in_trap)

        step_mag(in_trap) = ...
            abs(randn(sum(in_trap),1)) .* sigma_trap;

    end

    % Outside field: Cauchy movement resumes
    if any(outside)

        U = rand(sum(outside),1);

        step_mag(outside) = ...
            gamma .* tan((pi/2).*U);

    end

    %% Direction

    dir = zeros(Npart,1);

    % Attraction in cash crop

    x_clipped = ...
        max(0,min(zone_left,x_old));

    beta_space = ...
        (x_clipped ./ zone_left).^nu;

    beta_eff = ...
        beta .* beta_space;

    toward_prob = ...
        (1 + beta_eff)/2;

    r = rand(Npart,1);

    dir(was_in_cash) = -1;

    dir(was_in_cash & (r < toward_prob)) = 1;

    % Unbiased movement in trap crop and outside field

    unbiased_idx = ...
        in_trap | outside;

    if any(unbiased_idx)

        r2 = rand(sum(unbiased_idx),1);

        temp_dir = ones(sum(unbiased_idx),1);

        temp_dir(r2 > 0.5) = -1;

        dir(unbiased_idx) = temp_dir;

    end

    %% Position update

    x_new = x_old;

    x_new(active) = ...
        x_old(active) + ...
        dir(active).*step_mag(active);

    %% Entry-triggered removal

    crossed_into_trap = ...
        was_in_cash & ...
        (x_new > zone_left);

    to_kill_entry = ...
        crossed_into_trap & ...
        (rand(Npart,1) < alpha_entry);

    x = x_new;

    if any(to_kill_entry)

        is_killed(to_kill_entry) = true;
        x(to_kill_entry) = NaN;

    end

    %% Within-trap removal

    active = ~is_killed;

    in_strip = ...
        active & ...
        (x > zone_left) & ...
        (x <= L);

    to_kill_trap = ...
        in_strip & ...
        (rand(Npart,1) < alpha_step);

    if any(to_kill_trap)

        is_killed(to_kill_trap) = true;
        x(to_kill_trap) = NaN;

    end

    %% Save snapshots

    snap_idx = find(snap_steps == a);

    if ~isempty(snap_idx)

        x_snap_levy{snap_idx} = x;

    end

end

%% ==============================================================
%% SNAPSHOT PERCENTAGES
%% ==============================================================

n_snap = length(snap_days);

gauss_left    = zeros(n_snap,1);
gauss_cash    = zeros(n_snap,1);
gauss_trap    = zeros(n_snap,1);
gauss_right   = zeros(n_snap,1);
gauss_removed = zeros(n_snap,1);

levy_left     = zeros(n_snap,1);
levy_cash     = zeros(n_snap,1);
levy_trap     = zeros(n_snap,1);
levy_right    = zeros(n_snap,1);
levy_removed  = zeros(n_snap,1);

%% Short-distance Gaussian percentages

for j = 1:n_snap

    xj = x_snap_gauss{j};

    gauss_left(j) = ...
        100 * sum(xj < 0) / Npart;

    gauss_cash(j) = ...
        100 * sum((xj >= 0) & (xj <= zone_left)) / Npart;

    gauss_trap(j) = ...
        100 * sum((xj > zone_left) & (xj <= L)) / Npart;

    gauss_right(j) = ...
        100 * sum(xj > L) / Npart;

    gauss_removed(j) = ...
        100 * sum(isnan(xj)) / Npart;

end

%% Long-distance Cauchy percentages

for j = 1:n_snap

    xj = x_snap_levy{j};

    levy_left(j) = ...
        100 * sum(xj < 0) / Npart;

    levy_cash(j) = ...
        100 * sum((xj >= 0) & (xj <= zone_left)) / Npart;

    levy_trap(j) = ...
        100 * sum((xj > zone_left) & (xj <= L)) / Npart;

    levy_right(j) = ...
        100 * sum(xj > L) / Npart;

    levy_removed(j) = ...
        100 * sum(isnan(xj)) / Npart;

end

%% ==============================================================
%% PRINT SNAPSHOT PERCENTAGES
%% ==============================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf('SHORT-DISTANCE GAUSSIAN MOVEMENT\n');
fprintf('Percentages of initial population N = %d\n',Npart);
fprintf('===============================================================\n');
fprintf('Day    Left emig.   Cash crop   Trap crop   Right emig.   Removed   Total\n');
fprintf('--------------------------------------------------------------------------\n');

for j = 1:n_snap

    total_pct = ...
        gauss_left(j) + ...
        gauss_cash(j) + ...
        gauss_trap(j) + ...
        gauss_right(j) + ...
        gauss_removed(j);

    fprintf('%3.0f    %8.2f     %8.2f    %8.2f     %8.2f    %8.2f   %6.2f\n', ...
        snap_days(j), ...
        gauss_left(j), ...
        gauss_cash(j), ...
        gauss_trap(j), ...
        gauss_right(j), ...
        gauss_removed(j), ...
        total_pct);

end

fprintf('\n');
fprintf('===============================================================\n');
fprintf('LONG-DISTANCE CAUCHY MOVEMENT\n');
fprintf('Percentages of initial population N = %d\n',Npart);
fprintf('===============================================================\n');
fprintf('Day    Left emig.   Cash crop   Trap crop   Right emig.   Removed   Total\n');
fprintf('--------------------------------------------------------------------------\n');

for j = 1:n_snap

    total_pct = ...
        levy_left(j) + ...
        levy_cash(j) + ...
        levy_trap(j) + ...
        levy_right(j) + ...
        levy_removed(j);

    fprintf('%3.0f    %8.2f     %8.2f    %8.2f     %8.2f    %8.2f   %6.2f\n', ...
        snap_days(j), ...
        levy_left(j), ...
        levy_cash(j), ...
        levy_trap(j), ...
        levy_right(j), ...
        levy_removed(j), ...
        total_pct);

end

fprintf('\n');

%% ==============================================================
%% PLOT
%% ==============================================================

fig = figure('Color','w', ...
       'Position',[80 80 1200 650]);

%% Plot settings

bin_edges = -20:2.5:120;

bin_centres = ...
    bin_edges(1:end-1) + diff(bin_edges)/2;

cash_col   = [0.38 0.68 0.38];
trap_col   = [0.70 0.48 0.80];
emig_col   = [0.46 0.66 0.84];

trap_shade = [0.86 0.78 0.88];

%% Determine common percentage y-axis

max_y = 0;

for row = 1:2

    if row == 1
        snapshots = x_snap_gauss;
    else
        snapshots = x_snap_levy;
    end

    for k = 1:3

        data = snapshots{k};
        data = data(~isnan(data));

        counts = histcounts(data,bin_edges);

        perc_counts = ...
            100 * counts / Npart;

        max_y = max(max_y,max(perc_counts));

    end

end

common_ylim = [0 1.05*max_y];

%% Plot panels

for row = 1:2

    if row == 1
        snapshots = x_snap_gauss;
    else
        snapshots = x_snap_levy;
    end

    for k = 1:3

        panel_idx = ...
            (row-1)*3 + k;

        subplot(2,3,panel_idx);

        data = snapshots{k};
        data = data(~isnan(data));

        cash_data = ...
            data(data >= 0 & ...
                 data <= zone_left);

        trap_data = ...
            data(data > zone_left & ...
                 data <= L);

        emig_data = ...
            data(data < 0 | ...
                 data > L);

        %% Percentage in each spatial bin

        cash_counts = ...
            100 * histcounts(cash_data,bin_edges) / Npart;

        trap_counts = ...
            100 * histcounts(trap_data,bin_edges) / Npart;

        emig_counts = ...
            100 * histcounts(emig_data,bin_edges) / Npart;

        %% Trap-crop shading

        patch( ...
            [zone_left L L zone_left], ...
            [common_ylim(1) common_ylim(1) ...
             common_ylim(2) common_ylim(2)], ...
            trap_shade, ...
            'EdgeColor','none', ...
            'FaceAlpha',0.25);

        hold on;

        %% Bars

        bar( ...
            bin_centres,cash_counts,1, ...
            'FaceColor',cash_col, ...
            'EdgeColor','none');

        bar( ...
            bin_centres,trap_counts,1, ...
            'FaceColor',trap_col, ...
            'EdgeColor','none');

        bar( ...
            bin_centres,emig_counts,1, ...
            'FaceColor',emig_col, ...
            'EdgeColor','none');

        %% Boundaries

        plot( ...
            [0 0],common_ylim, ...
            'k--','LineWidth',1);

        plot( ...
            [zone_left zone_left], ...
            common_ylim, ...
            'k--','LineWidth',1);

        plot( ...
            [L L],common_ylim, ...
            'k--','LineWidth',1);

        %% Axes

        xlim([-20 120]);
        ylim(common_ylim);

        set(gca, ...
            'FontName','Arial', ...
            'FontSize',14, ...
            'LineWidth',1.0, ...
            'TickDir','out', ...
            'Box','on');

        %% Column titles

        if row == 1

            title( ...
                ['{\itt} = ' ...
                 num2str(snap_days(k)) ...
                 ' days'], ...
                'FontName','Arial', ...
                'FontSize',16, ...
                'FontWeight','normal');

        end

        %% x-axis label

        if row == 2

            xlabel( ...
                'space, {\itx} (m)', ...
                'FontName','Arial', ...
                'FontSize',16);

        end

        %% Row labels / y-axis labels

        if k == 1

            if row == 1

                ylabel( ...
                    {'(a) short-distance movement', ...
                     '% of initial population'}, ...
                    'FontName','Arial', ...
                    'FontSize',16);

            else

                ylabel( ...
                    {'(b) long-distance movement', ...
                     '% of initial population'}, ...
                    'FontName','Arial', ...
                    'FontSize',16);

            end

        end

        hold off;

    end

end

%% Export the PNG with a small white margin around all plotted content
if save_figure
    png_file = fullfile(fileparts(mfilename('fullpath')), 'figure_1.png');
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
