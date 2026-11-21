function figure_3()
% Removal after 20 days across trap allocation and attraction
% strength. Each panel subtracts the Gaussian cumulative removal percentage
% from a linked Cauchy case. Red indicates greater removal for the Cauchy
% case, blue indicates lower removal, and the contour marks estimated parity.
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
close all;

rng(1);  % Reproducible random sequence within each figure.

save_figure = true;


%% ==========================================================
% BASELINE PARAMETERS
% ===========================================================

dt = 0.5;                 % time step (h)

Npart = 10^3;             % initial population

num_days = 20;

p_steps = ...
    round(num_days*24/dt) + 1;

sigma = 0.5;              % Gaussian scale (m)

L = 100;                  % total field length (m)

%% ==========================================================
% ARRESTMENT IN TRAP CROP
% ===========================================================

rho = 1/3;

sigma_trap = ...
    rho*sigma;

%% ==========================================================
% REMOVAL
% ===========================================================

gamma_entry = 0.1;

gamma_day = 0.2;

gamma_dt = ...
    1 - (1-gamma_day).^(dt/24);

%% ==========================================================
% SPATIAL ATTRACTION PROFILE
% ===========================================================

nu = 1;

%% ==========================================================
% PARAMETER SWEEPS
% ===========================================================

% Use the same finer mesh on both axes.

mesh_step = 0.01;

% Trap-crop fraction:
%
% 0    = no trap crop
% 0.30 = 30% trap crop

trap_crop_frac_vals = ...
    0.01:mesh_step:0.30;

% Corresponding cash-crop fractions

k_vals = ...
    1-trap_crop_frac_vals;

% Maximum attraction strength

beta_vals = ...
    0.01:mesh_step:0.50;

%% ==========================================================
% LONG-DISTANCE RELATIONSHIPS
% ===========================================================

p_values = ...
    [0.1 0.5 0.721 0.9];

nP = ...
    length(p_values);

gamma_values = ...
    zeros(1,nP);

for ip = 1:nP

    p = ...
        p_values(ip);

    gamma_values(ip) = ...
        sqrt(2)*erfcinv(p) ./ ...
        tan((pi/2)*(1-p)) * sigma;

end

%% ==========================================================
% DISPLAY PARAMETERS
% ===========================================================

fprintf('\n');

fprintf('Gaussian sigma = %.3f m\n', ...
    sigma);

fprintf('Trap sigma = %.3f m\n', ...
    sigma_trap);

fprintf('Common parameter mesh = %.3f\n', ...
    mesh_step);

fprintf('\n');

for ip = 1:nP

    fprintf( ...
        'p = %.3f: gamma = %.3f m, gamma/sigma = %.3f\n', ...
        p_values(ip), ...
        gamma_values(ip), ...
        gamma_values(ip)/sigma);

end

fprintf('\n');

%% ==========================================================
% STORAGE
% ===========================================================

nB = ...
    length(beta_vals);

nK = ...
    length(k_vals);

% Gaussian removal:
%
% dimensions = trap fraction x beta

removed_gaussian = ...
    zeros(nK,nB);

% Long-distance removal:
%
% dimensions = trap fraction x beta x p

removed_long = ...
    zeros(nK,nB,nP);

% Difference:
%
% long-distance minus Gaussian

delta_removed = ...
    zeros(nK,nB,nP);

%% ==========================================================
% MAIN PARAMETER SWEEP
% ===========================================================

fprintf('Running trap-cropping heat-map simulations...\n\n');

for ik = 1:nK

    k = ...
        k_vals(ik);

    zone_left = ...
        k*L;

    fprintf( ...
        'Trap-crop fraction = %.3f\n', ...
        trap_crop_frac_vals(ik));

    for ib = 1:nB

        beta = ...
            beta_vals(ib);

        %% --------------------------------------------------
        % NO-TRAP BOUNDARY CONDITION
        % --------------------------------------------------
        %
        % When trap-crop fraction = 0, there is no trap crop
        % and hence no trap-associated removal under either
        % movement process.
        %
        % Therefore:
        %
        % R_short = R_long = Delta R = 0.
        %

        if trap_crop_frac_vals(ik) == 0

            removed_gaussian(ik,ib) = ...
                0;

            for ip = 1:nP

                removed_long(ik,ib,ip) = ...
                    0;

                delta_removed(ik,ib,ip) = ...
                    0;

            end

            continue;

        end

        %% --------------------------------------------------
        % COMMON INITIAL CONDITION
        % --------------------------------------------------

        x0 = ...
            zone_left*rand(Npart,1);

        %% --------------------------------------------------
        % SHORT-DISTANCE GAUSSIAN
        % --------------------------------------------------

        removed_gaussian(ik,ib) = ...
            simulate_trap_crop_removal( ...
            x0, ...
            false, ...
            NaN, ...
            sigma, ...
            sigma_trap, ...
            L, ...
            zone_left, ...
            beta, ...
            nu, ...
            gamma_entry, ...
            gamma_dt, ...
            Npart, ...
            p_steps);

        %% --------------------------------------------------
        % LONG-DISTANCE CASES
        % --------------------------------------------------

        for ip = 1:nP

            gamma = ...
                gamma_values(ip);

            removed_long(ik,ib,ip) = ...
                simulate_trap_crop_removal( ...
                x0, ...
                true, ...
                gamma, ...
                sigma, ...
                sigma_trap, ...
                L, ...
                zone_left, ...
                beta, ...
                nu, ...
                gamma_entry, ...
                gamma_dt, ...
                Npart, ...
                p_steps);

            delta_removed(ik,ib,ip) = ...
                removed_long(ik,ib,ip) - ...
                removed_gaussian(ik,ib);

        end

    end

end

%% ==========================================================
% REPORT RANGE
% ===========================================================

fprintf('\n');

for ip = 1:nP

    this_delta = ...
        delta_removed(:,:,ip);

    fprintf( ...
        ['p = %.3f: Delta removal range = %.2f to %.2f ' ...
         'percentage points\n'], ...
        p_values(ip), ...
        min(this_delta(:)), ...
        max(this_delta(:)));

end

fprintf('\n');

%% ==========================================================
% COMMON SYMMETRIC COLOUR SCALE
% ===========================================================

max_abs_delta = ...
    max(abs(delta_removed(:)));

% Round upward to a whole percentage point.

colour_limit = ...
    ceil(max_abs_delta);

if colour_limit < 1

    colour_limit = 1;

end

%% ==========================================================
% DIVERGING BLUE-WHITE-RED COLOUR MAP
% ===========================================================

nColours = 256;

halfN = ...
    nColours/2;

% Blue -> white

cmap1 = ...
    [ ...
    linspace(0,1,halfN)' ...
    linspace(0,1,halfN)' ...
    ones(halfN,1)];

% White -> red

cmap2 = ...
    [ ...
    ones(halfN,1) ...
    linspace(1,0,halfN)' ...
    linspace(1,0,halfN)'];

cmap = ...
    [cmap1; cmap2];

%% ==========================================================
% FIGURE
% ===========================================================

fig = figure( ...
    'Color','w', ...
    'Position',[70 60 1100 850]);

colormap(cmap);

panel_labels = ...
    {'(a)','(b)','(c)','(d)'};

ax = ...
    gobjects(nP,1);

for ip = 1:nP

    ax(ip) = ...
        subplot(2,2,ip);

    data_plot = ...
        delta_removed(:,:,ip);

    % Transpose data because:
    %
    % rows    = trap-crop fraction
    % columns = beta_max
    %
    % but the plot uses:
    %
    % x = trap-crop fraction
    % y = beta_max

    imagesc( ...
        trap_crop_frac_vals, ...
        beta_vals, ...
        data_plot.');

    set(gca, ...
        'YDir','normal');

    hold on;

    %% ------------------------------------------------------
    % ZERO-DIFFERENCE CONTOUR
    % ------------------------------------------------------

    if min(data_plot(:)) <= 0 && ...
       max(data_plot(:)) >= 0

        contour( ...
            trap_crop_frac_vals, ...
            beta_vals, ...
            data_plot.', ...
            [0 0], ...
            'k-', ...
            'LineWidth',1.5);

    end

    hold off;

    caxis( ...
        [-colour_limit colour_limit]);

    %% ------------------------------------------------------
    % AXIS LABELS
    % ------------------------------------------------------

    if ip == 3 || ip == 4

        xlabel( ...
            'trap-crop allocation (%)', ...
            'Interpreter','tex', ...
            'FontName','Arial', ...
            'FontSize',16);

    else

        xlabel('');

    end

    if ip == 1 || ip == 3

        ylabel( ...
            'maximum attraction strength, \beta_{max}', ...
            'Interpreter','tex', ...
            'FontName','Arial', ...
            'FontSize',16);

    else

        ylabel('');

    end

    %% ------------------------------------------------------
    % PANEL TITLE
    % ------------------------------------------------------

    if ip == 3

        title_text = ...
            sprintf( ...
            '%s p^{*} = %.3f', ...
            panel_labels{ip}, ...
            p_values(ip));

    else

        title_text = ...
            sprintf( ...
            '%s p = %.1f', ...
            panel_labels{ip}, ...
            p_values(ip));

    end

    title( ...
        title_text, ...
        'Interpreter','tex', ...
        'FontName','Arial', ...
        'FontSize',16, ...
        'FontWeight','normal');

    %% ------------------------------------------------------
    % AXIS FORMAT
    % ------------------------------------------------------

   set(gca, ...
    'XTick',[0.01 0.05 0.10 0.15 0.20 0.25 0.30], ...
    'XTickLabel',{'1','5','10','15','20','25','30'}, ...
        'YTick',[0.01 0.1 0.2 0.3 0.4 0.5], ...
        'FontName','Arial', ...
        'FontSize',14, ...
        'LineWidth',1.0, ...
        'TickDir','out', ...
        'Box','on', ...
        'XMinorTick','off', ...
        'YMinorTick','off');

    xlim( ...
        [0.01 0.30]);

    ylim( ...
        [0.01 0.50]);

    axis square;

end

%% ==========================================================
% SHARED COLOUR BAR
% ===========================================================

% Store subplot positions before adding colour bar.

ax_positions = ...
    cell(nP,1);

for ip = 1:nP

    ax_positions{ip} = ...
        get(ax(ip),'Position');

end

hcb = ...
    colorbar( ...
    'peer',ax(4));

set( ...
    hcb, ...
    'Position',[0.92 0.15 0.02 0.70], ...
    'FontName','Arial', ...
    'FontSize',14);

ylabel( ...
    hcb, ...
    '\Delta removal (% points)', ...
    'Interpreter','tex', ...
    'FontName','Arial', ...
    'FontSize',16);

% Restore subplot sizes.

for ip = 1:nP

    set( ...
        ax(ip), ...
        'Position',ax_positions{ip});

end

%% Export the PNG with a small white margin around all plotted content
if save_figure
    png_file = fullfile(fileparts(mfilename('fullpath')), 'figure_3.png');
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

fprintf('Simulation complete.\n');

end

%% =====================================================================
% SIMULATE ONE TRAP-CROPPING SCENARIO
% =====================================================================

function removed_pct = ...
    simulate_trap_crop_removal( ...
    x0, ...
    useLevy, ...
    gamma, ...
    sigma, ...
    sigma_trap, ...
    L, ...
    zone_left, ...
    beta, ...
    nu, ...
    gamma_entry, ...
    gamma_dt, ...
    Npart, ...
    p_steps)

%% ----------------------------------------------------------
% INITIAL CONDITION
% ----------------------------------------------------------

x = ...
    x0;

is_killed = ...
    false(Npart,1);

%% ----------------------------------------------------------
% DOES A TRAP CROP EXIST?
% ----------------------------------------------------------

trap_exists = ...
    zone_left < L;

%% ----------------------------------------------------------
% TIME LOOP
% ----------------------------------------------------------

for a = 2:p_steps

    active = ...
        ~is_killed;

    x_old = ...
        x;

    %% ------------------------------------------------------
    % LOCATION BEFORE MOVEMENT
    % ------------------------------------------------------

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
        ((x_old < 0) | ...
         (x_old > L));

    %% ------------------------------------------------------
    % STEP MAGNITUDES
    % ------------------------------------------------------

    step_mag = ...
        zeros(Npart,1);

    % Trap crop:
    %
    % BOTH dispersal processes use the same arrested
    % Gaussian movement.

    if any(in_trap)

        step_mag(in_trap) = ...
            abs(randn(sum(in_trap),1)) * ...
            sigma_trap;

    end

    % Cash crop and outside field:
    %
    % Gaussian or folded-Cauchy according to movement type.

    full_movement = ...
        was_in_cash | outside;

    if any(full_movement)

        nMove = ...
            sum(full_movement);

        if useLevy

            U = ...
                rand(nMove,1);

            step_mag(full_movement) = ...
                gamma * ...
                tan((pi/2)*U);

        else

            step_mag(full_movement) = ...
                abs(randn(nMove,1)) * ...
                sigma;

        end

    end

    %% ------------------------------------------------------
    % DIRECTION
    % ------------------------------------------------------

    dir_new = ...
        zeros(Npart,1);

    % Cash crop:
    %
    % Direction determined only by spatial attraction.
    % No directional persistence.

    if any(was_in_cash)

        x_cash = ...
            x_old(was_in_cash);

        beta_space = ...
            (x_cash./zone_left).^nu;

        beta_eff = ...
            beta*beta_space;

        toward_prob = ...
            (1+beta_eff)/2;

        toward = ...
            rand(sum(was_in_cash),1) < ...
            toward_prob;

        d = ...
            -ones(sum(was_in_cash),1);

        d(toward) = ...
            1;

        dir_new(was_in_cash) = ...
            d;

    end

    % Trap crop and outside:
    %
    % unbiased direction.

    unbiased = ...
        in_trap | outside;

    if any(unbiased)

        d = ...
            ones(sum(unbiased),1);

        d(rand(sum(unbiased),1) > 0.5) = ...
            -1;

        dir_new(unbiased) = ...
            d;

    end

    %% ------------------------------------------------------
    % POSITION UPDATE
    % ------------------------------------------------------

    x_new = ...
        x;

    x_new(active) = ...
        x_old(active) + ...
        dir_new(active).* ...
        step_mag(active);

    %% ------------------------------------------------------
    % ENTRY-TRIGGERED REMOVAL
    % ------------------------------------------------------

    % An individual starting in the cash crop is regarded as
    % having encountered/crossed the trap-crop boundary whenever
    % its new position lies beyond zone_left.
    %
    % This includes a long folded-Cauchy step that passes through
    % the entire trap-crop strip and lands beyond x = L.
    %
    % If trap crop fraction = 0, trap_exists = false and no
    % entry-triggered removal can occur.

    crossed_trap_boundary = ...
        trap_exists & ...
        was_in_cash & ...
        (x_new > zone_left);

    to_kill_entry = ...
        crossed_trap_boundary & ...
        (rand(Npart,1) < gamma_entry);

    if any(to_kill_entry)

        is_killed(to_kill_entry) = ...
            true;

        x_new(to_kill_entry) = ...
            NaN;

    end

    %% ------------------------------------------------------
    % WITHIN-TRAP REMOVAL
    % ------------------------------------------------------

    active_after_entry = ...
        ~is_killed;

    in_strip_after_move = ...
        trap_exists & ...
        active_after_entry & ...
        (x_new > zone_left) & ...
        (x_new <= L);

    to_kill_trap = ...
        in_strip_after_move & ...
        (rand(Npart,1) < gamma_dt);

    if any(to_kill_trap)

        is_killed(to_kill_trap) = ...
            true;

        x_new(to_kill_trap) = ...
            NaN;

    end

    %% ------------------------------------------------------
    % UPDATE POSITION
    % ------------------------------------------------------

    x = ...
        x_new;

end

%% ----------------------------------------------------------
% FINAL REMOVAL
% ----------------------------------------------------------

removed_pct = ...
    100*sum(is_killed)/Npart;

end
