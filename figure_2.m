clc; clear; close all;
rng(1);
%% ==========================================================
% PARAMETERS
% ===========================================================

dt      = 0.5;          % time step (h)
Npart   = 10^4;         % initial population
num_days = 20;

p_steps = round(num_days*24/dt) + 1;

sigma = 0.5;            % Gaussian scale (m)

L = 100;                % total field length (m)

trap_frac = 0.80;       % cash-crop fraction
zone_left = trap_frac * L;

%% ==========================================================
% ATTRACTION
% ===========================================================

beta = 0.1;             % maximum attraction strength
nu   = 1;               % spatial attraction profile

%% ==========================================================
% ARRESTMENT IN TRAP CROP
% ===========================================================

rho = 1/3;

sigma_trap = rho * sigma;

%% ==========================================================
% REMOVAL
% ===========================================================

gamma_entry = 0.1;

gamma_day = 0.2;

gamma_dt = ...
    1 - (1 - gamma_day).^(dt/24);

%% ==========================================================
% TIME
% ===========================================================

time_hours = (0:p_steps-1)*dt;
time_days  = time_hours/24;

%% ==========================================================
% p-VALUES AND DERIVED GAMMA
% ===========================================================

p1 = 0.1;

gamma1 = ...
    sqrt(2) * erfcinv(p1) ./ ...
    tan((pi/2)*(1-p1)) * sigma;

p2 = 0.5;

gamma2 = ...
    sqrt(2) * erfcinv(p2) ./ ...
    tan((pi/2)*(1-p2)) * sigma;

p3 = 0.721;             % optimal survival probability

gamma3 = ...
    sqrt(2) * erfcinv(p3) ./ ...
    tan((pi/2)*(1-p3)) * sigma;

p4 = 0.9;

gamma4 = ...
    sqrt(2) * erfcinv(p4) ./ ...
    tan((pi/2)*(1-p4)) * sigma;

gamma_values = ...
    [gamma1 gamma2 gamma3 gamma4];

fprintf('sigma = %.3f m\n',sigma);
fprintf('sigma_trap = %.3f m\n',sigma_trap);

fprintf('gamma (p = 0.1)   = %.3f m\n',gamma1);
fprintf('gamma (p = 0.5)   = %.3f m\n',gamma2);
fprintf('gamma (p*= 0.721) = %.3f m\n',gamma3);
fprintf('gamma (p = 0.9)   = %.3f m\n',gamma4);

%% ==========================================================
% STORAGE
% ===========================================================

nCases = 5;

cash_all = ...
    zeros(nCases,p_steps);

trap_all = ...
    zeros(nCases,p_steps);

removed_all = ...
    zeros(nCases,p_steps);

emigrated_all = ...
    zeros(nCases,p_steps);

%% ==========================================================
% COMMON INITIAL CONDITION
% ===========================================================

x0 = ...
    zone_left * rand(Npart,1);

%% ==========================================================
% SIMULATIONS
% ===========================================================

for icase = 1:nCases

    %% ------------------------------------------------------
    % Initial population
    % -------------------------------------------------------

    x = x0;

    is_killed = ...
        false(Npart,1);

    %% Initial population states

    active = ~is_killed;

    cash_all(icase,1) = ...
        100 * ...
        sum(active & ...
        (x >= 0) & ...
        (x <= zone_left)) / Npart;

    trap_all(icase,1) = ...
        100 * ...
        sum(active & ...
        (x > zone_left) & ...
        (x <= L)) / Npart;

    removed_all(icase,1) = ...
        100 * sum(is_killed) / Npart;

    emigrated_all(icase,1) = ...
        100 * ...
        sum(active & ...
        ((x < 0) | ...
         (x > L))) / Npart;

    %% ------------------------------------------------------
    % Time loop
    % -------------------------------------------------------

    for a = 2:p_steps

        active = ...
            ~is_killed;

        x_old = x;

        %% --------------------------------------------------
        % Current locations
        % ---------------------------------------------------

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

        %% --------------------------------------------------
        % Step magnitudes
        % ---------------------------------------------------

        step_mag = ...
            zeros(Npart,1);

        % --------------------------------
        % Cash crop
        % --------------------------------

        if icase == 1

            % Short-distance Gaussian movement

            if any(was_in_cash)

                step_mag(was_in_cash) = ...
                    abs(randn(sum(was_in_cash),1)) ...
                    .* sigma;

            end

        else

            % Long-distance Levy movement

            gamma = ...
                gamma_values(icase-1);

            if any(was_in_cash)

                U = ...
                    rand(sum(was_in_cash),1);

                step_mag(was_in_cash) = ...
                    gamma .* ...
                    tan((pi/2).*U);

            end

        end

        % --------------------------------
        % Trap crop
        % BOTH movement types switch to
        % the same arrested Gaussian walk
        % --------------------------------

        if any(in_trap)

            step_mag(in_trap) = ...
                abs(randn(sum(in_trap),1)) ...
                .* sigma_trap;

        end

        % --------------------------------
        % Outside field
        % --------------------------------

        if icase == 1

            % Gaussian movement resumes

            if any(outside)

                step_mag(outside) = ...
                    abs(randn(sum(outside),1)) ...
                    .* sigma;

            end

        else

            % Levy movement resumes

            gamma = ...
                gamma_values(icase-1);

            if any(outside)

                U = ...
                    rand(sum(outside),1);

                step_mag(outside) = ...
                    gamma .* ...
                    tan((pi/2).*U);

            end

        end

        %% --------------------------------------------------
        % Direction
        % ---------------------------------------------------

        dir = ...
            zeros(Npart,1);

        % --------------------------------
        % Attraction in cash crop
        % --------------------------------

        x_clipped = ...
            max(0,min(zone_left,x_old));

        beta_space = ...
            (x_clipped ./ zone_left).^nu;

        beta_eff = ...
            beta .* beta_space;

        toward_prob = ...
            (1 + beta_eff)/2;

        r = ...
            rand(Npart,1);

        dir(was_in_cash) = -1;

        dir(was_in_cash & ...
            (r < toward_prob)) = 1;

        % --------------------------------
        % Unbiased movement in trap crop
        % and outside field
        % --------------------------------

        unbiased_idx = ...
            in_trap | outside;

        if any(unbiased_idx)

            r2 = ...
                rand(sum(unbiased_idx),1);

            temp_dir = ...
                ones(sum(unbiased_idx),1);

            temp_dir(r2 > 0.5) = -1;

            dir(unbiased_idx) = ...
                temp_dir;

        end

        %% --------------------------------------------------
        % Position update
        % ---------------------------------------------------

        x_new = x_old;

        x_new(active) = ...
            x_old(active) + ...
            dir(active).*step_mag(active);

        %% --------------------------------------------------
        % Entry-triggered removal
        % ---------------------------------------------------

        crossed_into_trap = ...
            was_in_cash & ...
            (x_new > zone_left);

        to_kill_entry = ...
            crossed_into_trap & ...
            (rand(Npart,1) < gamma_entry);

        x = x_new;

        if any(to_kill_entry)

            is_killed(to_kill_entry) = true;

            x(to_kill_entry) = NaN;

        end

        %% --------------------------------------------------
        % Within-trap removal
        % ---------------------------------------------------

        active = ...
            ~is_killed;

        in_strip = ...
            active & ...
            (x > zone_left) & ...
            (x <= L);

        to_kill_trap = ...
            in_strip & ...
            (rand(Npart,1) < gamma_dt);

        if any(to_kill_trap)

            is_killed(to_kill_trap) = true;

            x(to_kill_trap) = NaN;

        end

        %% --------------------------------------------------
        % Population states
        % ---------------------------------------------------

        active = ...
            ~is_killed;

        cash_all(icase,a) = ...
            100 * ...
            sum(active & ...
            (x >= 0) & ...
            (x <= zone_left)) / Npart;

        trap_all(icase,a) = ...
            100 * ...
            sum(active & ...
            (x > zone_left) & ...
            (x <= L)) / Npart;

        removed_all(icase,a) = ...
            100 * ...
            sum(is_killed) / Npart;

        emigrated_all(icase,a) = ...
            100 * ...
            sum(active & ...
            ((x < 0) | ...
             (x > L))) / Npart;

    end

end

%% ==========================================================
% PRINT RESULTS FOR FIGURE 2
% ===========================================================

% Find indices corresponding to days 10 and 20

[~,idx10] = min(abs(time_days - 10));
[~,idx20] = min(abs(time_days - 20));

case_names = { ...
    'Short-distance Gaussian'; ...
    'Long-distance p = 0.1'; ...
    'Long-distance p = 0.5'; ...
    'Long-distance p* = 0.721'; ...
    'Long-distance p = 0.9'};

fprintf('\n');
fprintf('============================================================\n');
fprintf('FIGURE 2 RESULTS SUMMARY\n');
fprintf('All values are %% of the initial population N = %d\n',Npart);
fprintf('============================================================\n');

for icase = 1:nCases

    [peak_trap,idx_peak_trap] = ...
        max(trap_all(icase,:));

    peak_trap_day = ...
        time_days(idx_peak_trap);

    [peak_emig,idx_peak_emig] = ...
        max(emigrated_all(icase,:));

    peak_emig_day = ...
        time_days(idx_peak_emig);

    fprintf('\n');
    fprintf('%s\n',case_names{icase});
    fprintf('------------------------------------------------------------\n');

    fprintf('Day 10:\n');
    fprintf('  Cash crop   = %6.2f%%\n',cash_all(icase,idx10));
    fprintf('  Trap crop   = %6.2f%%\n',trap_all(icase,idx10));
    fprintf('  Removed     = %6.2f%%\n',removed_all(icase,idx10));
    fprintf('  Emigrated   = %6.2f%%\n',emigrated_all(icase,idx10));

    fprintf('Day 20:\n');
    fprintf('  Cash crop   = %6.2f%%\n',cash_all(icase,idx20));
    fprintf('  Trap crop   = %6.2f%%\n',trap_all(icase,idx20));
    fprintf('  Removed     = %6.2f%%\n',removed_all(icase,idx20));
    fprintf('  Emigrated   = %6.2f%%\n',emigrated_all(icase,idx20));

    fprintf('Peak trap-crop occupancy = %6.2f%% at day %.2f\n', ...
        peak_trap,peak_trap_day);

    fprintf('Peak emigration          = %6.2f%% at day %.2f\n', ...
        peak_emig,peak_emig_day);

end

%% ==========================================================
% COMPACT COMPARISON TABLE
% ===========================================================

fprintf('\n');
fprintf('==============================================================================================\n');
fprintf('COMPACT COMPARISON\n');
fprintf('==============================================================================================\n');
fprintf('Case                         Cash20   Trap20   Removed20   Emig20   PeakTrap   PeakTrapDay\n');
fprintf('----------------------------------------------------------------------------------------------\n');

for icase = 1:nCases

    [peak_trap,idx_peak_trap] = ...
        max(trap_all(icase,:));

    peak_trap_day = ...
        time_days(idx_peak_trap);

    fprintf('%-28s %7.2f   %6.2f     %7.2f   %6.2f     %7.2f       %6.2f\n', ...
        case_names{icase}, ...
        cash_all(icase,idx20), ...
        trap_all(icase,idx20), ...
        removed_all(icase,idx20), ...
        emigrated_all(icase,idx20), ...
        peak_trap, ...
        peak_trap_day);

end

fprintf('==============================================================================================\n');
fprintf('\n');

%% ==========================================================
% PLOT
% ===========================================================

figure('Color','w', ...
       'Position',[100 100 1000 800]);

%% ==========================================================
% PANEL (a): CASH-CROP OCCUPANCY
% ===========================================================

subplot(2,2,1);
hold on;

plot(time_days, ...
    cash_all(1,:), '-', ...
    'Color','k', ...
    'LineWidth',1.8);

plot(time_days, ...
    cash_all(2,:), '--', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(time_days, ...
    cash_all(3,:), '-.', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(time_days, ...
    cash_all(4,:), '-', ...
    'Color','r', ...
    'LineWidth',1.8);

plot(time_days, ...
    cash_all(5,:), ':', ...
    'Color','b', ...
    'LineWidth',1.8);

ylabel('cash-crop occupancy (%)', ...
    'FontSize',16, ...
    'FontName','Arial');

title('(a)', ...
    'FontSize',16, ...
    'FontName','Arial', ...
    'FontWeight','normal');

box on;

xlim([0 num_days]);
ylim([0 100]);

set(gca,'XTick',[0 5 10 15 20 25 30]);

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out');

hold off;

%% ==========================================================
% PANEL (b): TRAP-CROP OCCUPANCY
% ===========================================================

subplot(2,2,2);
hold on;

h1 = plot(time_days, ...
    trap_all(1,:), '-', ...
    'Color','k', ...
    'LineWidth',1.8);

h2 = plot(time_days, ...
    trap_all(2,:), '--', ...
    'Color','b', ...
    'LineWidth',1.8);

h3 = plot(time_days, ...
    trap_all(3,:), '-.', ...
    'Color','b', ...
    'LineWidth',1.8);

h4 = plot(time_days, ...
    trap_all(4,:), '-', ...
    'Color','r', ...
    'LineWidth',1.8);

h5 = plot(time_days, ...
    trap_all(5,:), ':', ...
    'Color','b', ...
    'LineWidth',1.8);

% Dummy invisible handles used as legend headings

hHeaderShort = ...
    plot(NaN,NaN, ...
    'LineStyle','none', ...
    'Marker','none', ...
    'Color','w');

hHeaderLong = ...
    plot(NaN,NaN, ...
    'LineStyle','none', ...
    'Marker','none', ...
    'Color','w');

ylabel('trap-crop occupancy (%)', ...
    'FontSize',16, ...
    'FontName','Arial');

title('(b)', ...
    'FontSize',16, ...
    'FontName','Arial', ...
    'FontWeight','normal');

lgd = legend( ...
    [hHeaderShort h1 hHeaderLong h2 h3 h4 h5], { ...
    'short-distance dispersal', ...
    sprintf('   \\sigma = %.1f', sigma), ...
    'long-distance dispersal', ...
    sprintf('   \\color{blue}{\\itp} = %.1f', p1), ...
    sprintf('   \\color{blue}{\\itp} = %.1f', p2), ...
    sprintf('   \\color{red}{\\itp}^{*} = %.3f', p3), ...
    sprintf('   \\color{blue}{\\itp} = %.1f', p4)}, ...
   'Location','northwest', ...
    'Interpreter','tex', ...
    'FontSize',14, ...
    'FontName','Arial', ...
    'Box','off');

% Match legend text colours to trajectories

legend_text = findobj(lgd,'Type','text');

for ii = 1:length(legend_text)

    txt = get(legend_text(ii),'String');

    if ~isempty(strfind(txt,'0.721'))

        % Optimal case
        set(legend_text(ii),'Color','r');

    elseif ~isempty(strfind(txt,'long-distance')) || ...
           ~isempty(strfind(txt,'0.1')) || ...
           ~isempty(strfind(txt,'0.5')) || ...
           ~isempty(strfind(txt,'0.9'))

        % Other long-distance cases
        set(legend_text(ii),'Color','b');

    else

        % Gaussian / short-distance case
        set(legend_text(ii),'Color','k');

    end

end

box on;

xlim([0 num_days]);
ylim([0 100]);

set(gca,'XTick',[0 5 10 15 20 25 30]);

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out');

hold off;

%% ==========================================================
% PANEL (c): REMOVED
% ===========================================================

subplot(2,2,3);
hold on;

plot(time_days, ...
    removed_all(1,:), '-', ...
    'Color','k', ...
    'LineWidth',1.8);

plot(time_days, ...
    removed_all(2,:), '--', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(time_days, ...
    removed_all(3,:), '-.', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(time_days, ...
    removed_all(4,:), '-', ...
    'Color','r', ...
    'LineWidth',1.8);

plot(time_days, ...
    removed_all(5,:), ':', ...
    'Color','b', ...
    'LineWidth',1.8);

xlabel('time (days)', ...
    'FontSize',16, ...
    'FontName','Arial');

ylabel('removed (%)', ...
    'FontSize',16, ...
    'FontName','Arial');

title('(c)', ...
    'FontSize',16, ...
    'FontName','Arial', ...
    'FontWeight','normal');

box on;

xlim([0 num_days]);
ylim([0 100]);

set(gca,'XTick',[0 5 10 15 20 25 30]);

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out');

hold off;

%% ==========================================================
% PANEL (d): EMIGRATED
% ===========================================================

subplot(2,2,4);
hold on;

plot(time_days, ...
    emigrated_all(1,:), '-', ...
    'Color','k', ...
    'LineWidth',1.8);

plot(time_days, ...
    emigrated_all(2,:), '--', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(time_days, ...
    emigrated_all(3,:), '-.', ...
    'Color','b', ...
    'LineWidth',1.8);

plot(time_days, ...
    emigrated_all(4,:), '-', ...
    'Color','r', ...
    'LineWidth',1.8);

plot(time_days, ...
    emigrated_all(5,:), ':', ...
    'Color','b', ...
    'LineWidth',1.8);

xlabel('time (days)', ...
    'FontSize',16, ...
    'FontName','Arial');

ylabel('emigrated (%)', ...
    'FontSize',16, ...
    'FontName','Arial');

title('(d)', ...
    'FontSize',16, ...
    'FontName','Arial', ...
    'FontWeight','normal');

box on;

xlim([0 num_days]);
ylim([0 100]);

set(gca,'XTick',[0 5 10 15 20 25 30]);

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out');

hold off;