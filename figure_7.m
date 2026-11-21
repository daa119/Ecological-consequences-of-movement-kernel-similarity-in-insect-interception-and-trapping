function figure_7()
% An x-y projection of the three-dimensional trajectories in
% Figure 6. Separate runs reproduce the same draws with the same seed
% and parameters. In panel (c), long steps may leave and re-enter the
% displayed region. The plotted window does not constrain movement.
%
% The Chi(3) scale is calibrated from nominal flight speed.
% The folded-Cauchy model has no finite mean step length. Until the first
% crossing of the panel plane, headings are biased toward an individually
% assigned target on the square panel. Later ordinary headings are isotropic.
% Capture simulations test the whole segment for panel contact before
% checking endpoint escape. After failed adhesion, the trajectory resumes
% from the contact point with a small offset and an outward heading.
% In the free-flight panels, crossing the plane switches off target
% attraction and the trajectory continues without a contact response.

clear; clc; close all;
rng(1);

save_figure = true;

N = 10;

%% ==========================================================
% ECOLOGICAL PARAMETERS
% ===========================================================

dt = 0.02;                 % time step (s)

meanSpeed = 25;            % mean flight speed (cm/s)

meanStep = ...
    meanSpeed * dt;        % mean step length (cm)

% ----------------------------------------------------------
% For a 3D Gaussian random walk, radial step length follows
% a Chi(3) distribution:
%
% E[l] = 2*sigma*sqrt(2/pi)
%
% Hence:
% ----------------------------------------------------------

sigma = ...
    meanStep * sqrt(2*pi) / 4;

maxSteps = 1000;            % 20 s total duration

%% ==========================================================
% OPTIMAL 3D GAUSSIAN-CAUCHY RELATIONSHIP
% ===========================================================

p_opt = 0.650;

gamma_ratio = 2.089;

gamma = ...
    gamma_ratio * sigma;

%% ==========================================================
% DIRECTIONAL BIAS
% ===========================================================

kappa_b = 0.3;  % vMF concentration before the first plane crossing.

%% ==========================================================
% DISPLAY PARAMETER VALUES
% ===========================================================

fprintf('\n');

fprintf('Mean flight speed = %.3f cm/s\n', ...
    meanSpeed);

fprintf('Time step = %.3f s\n', ...
    dt);

fprintf('Mean step length = %.3f cm\n', ...
    meanStep);

fprintf('Gaussian scale sigma = %.3f cm\n', ...
    sigma);

fprintf('\n');

fprintf('Optimal p = %.3f\n', ...
    p_opt);

fprintf('gamma/sigma = %.3f\n', ...
    gamma_ratio);

fprintf('Cauchy scale gamma = %.3f cm\n', ...
    gamma);

fprintf('\n');

fprintf('Directional bias concentration kappa_b = %.3f\n', ...
    kappa_b);

fprintf('\n');

%% ==========================================================
% GEOMETRY
% ===========================================================

% Square interception panel

W = 10;                    % width (cm)
H = 10;                    % height (cm)

% Circular release region

xStart = -50;              % cm

releaseRad = 10;           % cm

%% ==========================================================
% PANEL CONTACT
% ===========================================================

pStick = 0.5;

kappa_post = 5.0;

epsShift = 1e-3;           % cm

%% ==========================================================
% DISPLAY / DOMAIN
% ===========================================================

nInitShow = N;

nTrajShow = N;

xEscapePos = 100;

xEscapeNeg = -100;

rEscape = 100;

xPlotMin = -60;

xPlotMax = 20;

yzPlot = 40;

% Trajectory colours for panels (b) and (c)
% Green is deliberately excluded because the release disc is green.

trajColors = ...
    [ ...
    1.00 0.55 0.00; ...   % orange
    0.60 0.30 0.80; ...   % purple
    0.20 0.60 0.90; ...   % blue
    0.90 0.40 0.55; ...   % pink
    0.95 0.75 0.20; ...   % yellow
    0.55 0.35 0.20; ...   % brown
    0.30 0.75 0.85; ...   % cyan
    0.85 0.25 0.25; ...   % red
    0.50 0.50 0.50; ...   % grey
    0.25 0.25 0.65];      % indigo

%% ==========================================================
% INITIAL RELEASE POINTS
% ===========================================================

pts = ...
    sample_release_points( ...
    nInitShow, ...
    xStart, ...
    releaseRad);

%% ==========================================================
% SHORT-DISTANCE TRAJECTORIES
% ===========================================================

paths_srw = ...
    cell(nTrajShow,1);

for i = 1:nTrajShow

    paths_srw{i} = ...
        simulate_one_trajectory_free( ...
        sigma, ...
        gamma, ...
        maxSteps, ...
        xStart, ...
        releaseRad, ...
        kappa_b, ...
        xEscapePos, ...
        xEscapeNeg, ...
        rEscape, ...
        W, ...
        H, ...
        false);

end

%% ==========================================================
% LONG-DISTANCE TRAJECTORIES
% ===========================================================

paths_lw = ...
    cell(nTrajShow,1);

for i = 1:nTrajShow

    paths_lw{i} = ...
        simulate_one_trajectory_free( ...
        sigma, ...
        gamma, ...
        maxSteps, ...
        xStart, ...
        releaseRad, ...
        kappa_b, ...
        xEscapePos, ...
        xEscapeNeg, ...
        rEscape, ...
        W, ...
        H, ...
        true);

end

%% ==========================================================
% REPRESENTATIVE GAUSSIAN PANEL INTERCEPTION
% ===========================================================

foundGauss = false;

nTry = 0;

while ~foundGauss && nTry < 50000

    nTry = nTry + 1;

    [status_gauss,path] = ...
        simulate_one_trajectory_with_status_contact( ...
        sigma, ...
        gamma, ...
        maxSteps, ...
        xStart, ...
        releaseRad, ...
        W, ...
        H, ...
        kappa_b, ...
        xEscapePos, ...
        xEscapeNeg, ...
        rEscape, ...
        pStick, ...
        kappa_post, ...
        epsShift, ...
        false);

    if status_gauss == 1

        path_gauss_intercept = path;

        foundGauss = true;

    end

end

%% ==========================================================
% REPRESENTATIVE LEVY PANEL INTERCEPTION
% ===========================================================

foundLevy = false;

nTry = 0;

while ~foundLevy && nTry < 50000

    nTry = nTry + 1;

    [status_levy,path] = ...
        simulate_one_trajectory_with_status_contact( ...
        sigma, ...
        gamma, ...
        maxSteps, ...
        xStart, ...
        releaseRad, ...
        W, ...
        H, ...
        kappa_b, ...
        xEscapePos, ...
        xEscapeNeg, ...
        rEscape, ...
        pStick, ...
        kappa_post, ...
        epsShift, ...
        true);

    if status_levy == 1

        path_levy_intercept = path;

        foundLevy = true;

    end

end

%% ==========================================================
% FIGURE 7: 2D x-y PROJECTIONS
% ===========================================================

fig7 = figure( ...
    'Color','w', ...
    'Position',[100 80 1000 800]);

%% ----------------------------------------------------------
% PANEL (a): INITIAL CONFIGURATION
% ----------------------------------------------------------

subplot(2,2,1);

hold on;

draw_panel_xy(W);

draw_release_segment_xy( ...
    xStart, ...
    releaseRad);

plot( ...
    pts(:,1), ...
    pts(:,2), ...
    'ko', ...
    'MarkerSize',5, ...
    'MarkerFaceColor',[0.2 0.2 0.2]);

title( ...
    '(a) initial configuration', ...
    'FontName','Arial', ...
    'FontSize',16, ...
    'FontWeight','normal');

format_axes_xy( ...
    xPlotMin, ...
    xPlotMax, ...
    yzPlot);

hold off;

%% ----------------------------------------------------------
% PANEL (b): SHORT-DISTANCE movement
% ----------------------------------------------------------

subplot(2,2,2);

hold on;

draw_panel_xy(W);

draw_release_segment_xy( ...
    xStart, ...
    releaseRad);

set(gca, ...
    'ColorOrder',trajColors, ...
    'ColorOrderIndex',1);

for i = 1:nTrajShow

    path = ...
        paths_srw{i};

    plot( ...
        path(:,1), ...
        path(:,2), ...
        '-', ...
        'LineWidth',1.2);

end

title( ...
    '(b) short-distance movement', ...
    'FontName','Arial', ...
    'FontSize',16, ...
    'FontWeight','normal');

format_axes_xy( ...
    xPlotMin, ...
    xPlotMax, ...
    yzPlot);

hold off;

%% ----------------------------------------------------------
% PANEL (c): LONG-DISTANCE movement
% ----------------------------------------------------------

subplot(2,2,3);

hold on;

draw_panel_xy(W);

draw_release_segment_xy( ...
    xStart, ...
    releaseRad);

set(gca, ...
    'ColorOrder',trajColors, ...
    'ColorOrderIndex',1);

for i = 1:nTrajShow

    path = ...
        paths_lw{i};

    plot( ...
        path(:,1), ...
        path(:,2), ...
        '-', ...
        'LineWidth',1.2);

end

title( ...
    '(c) long-distance movement', ...
    'FontName','Arial', ...
    'FontSize',16, ...
    'FontWeight','normal');

format_axes_xy( ...
    xPlotMin, ...
    xPlotMax, ...
    yzPlot);

hold off;

%% ----------------------------------------------------------
% PANEL (d): REPRESENTATIVE PANEL INTERCEPTION
% ----------------------------------------------------------

subplot(2,2,4);

hold on;

draw_panel_xy(W);

draw_release_segment_xy( ...
    xStart, ...
    releaseRad);

if foundGauss

    hGauss = ...
        plot( ...
        path_gauss_intercept(:,1), ...
        path_gauss_intercept(:,2), ...
        'k-', ...
        'LineWidth',1.2);

else

    hGauss = ...
        plot(NaN,NaN, ...
        'k-', ...
        'LineWidth',1.2);

end

if foundLevy

    hLevy = ...
        plot( ...
        path_levy_intercept(:,1), ...
        path_levy_intercept(:,2), ...
        'b-', ...
        'LineWidth',1.2);

else

    hLevy = ...
        plot(NaN,NaN, ...
        'b-', ...
        'LineWidth',1.2);

end

title( ...
    '(d) panel interception', ...
    'FontName','Arial', ...
    'FontSize',16, ...
    'FontWeight','normal');

legend( ...
    [hGauss hLevy], ...
    { ...
    'short-distance movement', ...
    'long-distance movement'}, ...
    'Location','southeast', ...
    'FontName','Arial', ...
    'FontSize',12, ...
    'Box','off');

format_axes_xy( ...
    xPlotMin, ...
    xPlotMax, ...
    yzPlot);

hold off;

%% Export the PNG with a small white margin around all plotted content
if save_figure
    png_file = fullfile(fileparts(mfilename('fullpath')), 'figure_7.png');
    dpi = 300;
    drawnow nocallbacks;
    if ~isgraphics(fig7, 'figure')
        warning('Figure 7 closed before export. figure_7.png was not updated.');
        return;
    end
    set(fig7, 'PaperPositionMode', 'auto');
    set(0, 'CurrentFigure', fig7);
    print(png_file, '-dpng', sprintf('-r%d',dpi));

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

end

%% =====================================================================
% FREE TRAJECTORY
% =====================================================================

function path = simulate_one_trajectory_free( ...
    sigma, ...
    gamma, ...
    maxSteps, ...
    xStart, ...
    releaseRad, ...
    kappa_b, ...
    xEscapePos, ...
    xEscapeNeg, ...
    rEscape, ...
    W, ...
    H, ...
    useLevy)

%% Initial position

theta = ...
    2*pi*rand;

rr = ...
    releaseRad*sqrt(rand);

y0 = ...
    rr*cos(theta);

z0 = ...
    rr*sin(theta);

r = ...
    [xStart y0 z0];

%% Random target on panel

yTar = ...
    -W/2 + W*rand;

zTar = ...
    -H/2 + H*rand;

pTar = ...
    [0 yTar zTar];

attractionOn = true;

%% Storage

path = ...
    zeros(maxSteps+1,3);

idx = 1;

path(idx,:) = r;

%% Movement loop

for t = 1:maxSteps

    %% ------------------------------------------------------
    % DIRECTION
    % -------------------------------------------------------

    if attractionOn

        g = ...
            pTar-r;

        ng = ...
            norm(g);

        if ng < 1e-14

            g = ...
                [1 0 0];

        else

            g = ...
                g/ng;

        end

        v_new = ...
            sample_vMF_3d( ...
            g, ...
            kappa_b);

    else

        v_new = ...
            rand_unit_vector();

    end

    %% ------------------------------------------------------
    % STEP LENGTH
    % -------------------------------------------------------

    if useLevy

        ell = ...
            sample_folded_cauchy( ...
            gamma);

    else

        ell = ...
            sample_chi3( ...
            sigma);

    end

    %% ------------------------------------------------------
    % POSITION UPDATE
    % -------------------------------------------------------

    r_new = ...
        r + ell*v_new;

    %% ------------------------------------------------------
    % REMOVE BIAS AFTER CROSSING PANEL PLANE
    % -------------------------------------------------------

    if attractionOn

        [~,crossedPlane,~] = ...
            segment_hits_panel( ...
            r, ...
            r_new, ...
            W, ...
            H);

        if crossedPlane

            attractionOn = false;

        end

    end

    idx = ...
        idx+1;

    path(idx,:) = ...
        r_new;

    %% ------------------------------------------------------
    % ESCAPE
    % -------------------------------------------------------

    if r_new(1) > xEscapePos || ...
       r_new(1) < xEscapeNeg || ...
       norm(r_new) > rEscape

        path = ...
            path(1:idx,:);

        return;

    end

    r = ...
        r_new;

end

path = ...
    path(1:idx,:);

end

%% =====================================================================
% TRAJECTORY WITH PANEL CONTACT
% =====================================================================

function [status,path] = ...
    simulate_one_trajectory_with_status_contact( ...
    sigma, ...
    gamma, ...
    maxSteps, ...
    xStart, ...
    releaseRad, ...
    W, ...
    H, ...
    kappa_b, ...
    xEscapePos, ...
    xEscapeNeg, ...
    rEscape, ...
    pStick, ...
    kappa_post, ...
    epsShift, ...
    useLevy)

% status:
%
%   1 = captured on panel
%   2 = crossed panel plane outside panel bounds
%   3 = no interception before termination

%% Initial position

theta = ...
    2*pi*rand;

rr = ...
    releaseRad*sqrt(rand);

y0 = ...
    rr*cos(theta);

z0 = ...
    rr*sin(theta);

r = ...
    [xStart y0 z0];

%% Random target on panel

yTar = ...
    -W/2 + W*rand;

zTar = ...
    -H/2 + H*rand;

pTar = ...
    [0 yTar zTar];

attractionOn = true;

%% Storage

path = ...
    zeros(3*maxSteps+20,3);

idx = 1;

path(idx,:) = ...
    r;

leaked = false;

%% Movement

for t = 1:maxSteps

    %% ------------------------------------------------------
    % DIRECTION
    % -------------------------------------------------------

    if attractionOn

        g = ...
            pTar-r;

        ng = ...
            norm(g);

        if ng < 1e-14

            g = ...
                [1 0 0];

        else

            g = ...
                g/ng;

        end

        v_new = ...
            sample_vMF_3d( ...
            g, ...
            kappa_b);

    else

        v_new = ...
            rand_unit_vector();

    end

    %% ------------------------------------------------------
    % STEP LENGTH
    % -------------------------------------------------------

    if useLevy

        ell = ...
            sample_folded_cauchy( ...
            gamma);

    else

        ell = ...
            sample_chi3( ...
            sigma);

    end

    r_new = ...
        r + ell*v_new;

    %% ------------------------------------------------------
    % PANEL INTERSECTION
    % -------------------------------------------------------

    [isHit,crossedPlane,cp] = ...
        segment_hits_panel( ...
        r, ...
        r_new, ...
        W, ...
        H);

    %% ------------------------------------------------------
    % PANEL CONTACT
    % -------------------------------------------------------

    if isHit

        idx = ...
            idx+1;

        path(idx,:) = ...
            cp;

        if rand < pStick

            status = 1;

            path = ...
                path(1:idx,:);

            return;

        else

            %% Failed adhesion

            v_post = ...
                sample_post_contact_heading( ...
                kappa_post);

            if useLevy

                ell_post = ...
                    sample_folded_cauchy( ...
                    gamma);

            else

                ell_post = ...
                    sample_chi3( ...
                    sigma);

            end

            r = ...
                cp + ...
                epsShift*v_post;

            r_new2 = ...
                r + ...
                ell_post*v_post;

            idx = ...
                idx+1;

            path(idx,:) = ...
                r_new2;

            if r_new2(1) > xEscapePos || ...
               r_new2(1) < xEscapeNeg || ...
               norm(r_new2) > rEscape

                status = 3;

                path = ...
                    path(1:idx,:);

                return;

            end

            r = ...
                r_new2;

            attractionOn = false;

            continue;

        end

    %% ------------------------------------------------------
    % CROSS PANEL PLANE OUTSIDE PANEL
    % -------------------------------------------------------

    elseif crossedPlane

        % Keep the endpoint of the original step. This plane crossing
        % switches off target attraction, with no contact or rebound.
        leaked = true;
        attractionOn = false;

        idx = idx + 1;
        path(idx,:) = r_new;

    else

        idx = ...
            idx+1;

        path(idx,:) = ...
            r_new;

    end

    %% ------------------------------------------------------
    % ESCAPE
    % -------------------------------------------------------

    if r_new(1) > xEscapePos || ...
       r_new(1) < xEscapeNeg || ...
       norm(r_new) > rEscape

        if leaked

            status = 2;

        else

            status = 3;

        end

        path = ...
            path(1:idx,:);

        return;

    end

    r = ...
        r_new;

end

%% Final status

if leaked

    status = 2;

else

    status = 3;

end

path = ...
    path(1:idx,:);

end

%% =====================================================================
% CHI(3) STEP LENGTH
% =====================================================================

function ell = sample_chi3(sigma)

ell = ...
    sigma * ...
    sqrt( ...
    randn^2 + ...
    randn^2 + ...
    randn^2);

end

%% =====================================================================
% FOLDED-CAUCHY STEP LENGTH
% =====================================================================

function ell = sample_folded_cauchy(gamma)

U = ...
    rand;

ell = ...
    gamma * ...
    tan((pi/2)*U);

end

%% =====================================================================
% POST-CONTACT DIRECTION
% =====================================================================

function v = ...
    sample_post_contact_heading( ...
    kappa_post)

mu_post = ...
    [-1 0 0];

while true

    v = ...
        sample_vMF_3d( ...
        mu_post, ...
        kappa_post);

    if v(1) < 0

        return;

    end

end

end

%% =====================================================================
% PANEL INTERSECTION
% =====================================================================

function [tf,crossedPlane,point] = ...
    segment_hits_panel( ...
    r1, ...
    r2, ...
    W, ...
    H)

x1 = ...
    r1(1);

x2 = ...
    r2(1);

crossedPlane = false;

point = ...
    [NaN NaN NaN];

if (x1 < 0 && x2 < 0) || ...
   (x1 > 0 && x2 > 0)

    tf = false;

    return;

end

dx = ...
    x2-x1;

if abs(dx) < 1e-14

    tf = false;

    return;

end

tau = ...
    -x1/dx;

if tau < 0 || tau > 1

    tf = false;

    return;

end

crossedPlane = true;

yint = ...
    r1(2) + ...
    tau*(r2(2)-r1(2));

zint = ...
    r1(3) + ...
    tau*(r2(3)-r1(3));

point = ...
    [0 yint zint];

tf = ...
    (abs(yint) <= W/2) && ...
    (abs(zint) <= H/2);

end

%% =====================================================================
% VON MISES-FISHER DIRECTION
% =====================================================================

function x = ...
    sample_vMF_3d( ...
    mu, ...
    kappa)

mu = ...
    mu(:)'/norm(mu);

if kappa < 1e-10

    x = ...
        rand_unit_vector();

    return;

end

u = ...
    rand;

a = ...
    exp(-kappa);

b = ...
    exp(kappa);

w = ...
    (1/kappa) * ...
    log(a + u*(b-a));

w = ...
    min(1,max(-1,w));

phi = ...
    2*pi*rand;

s = ...
    sqrt(max(0,1-w^2));

xp = ...
    [ ...
    s*cos(phi), ...
    s*sin(phi), ...
    w];

north = ...
    [0 0 1];

if norm(mu-north) < 1e-12

    x = ...
        xp;

elseif norm(mu+north) < 1e-12

    x = ...
        [xp(1) xp(2) -xp(3)];

else

    R = ...
        rotation_matrix_from_a_to_b( ...
        north, ...
        mu);

    x = ...
        (R*xp(:))';

    x = ...
        x/norm(x);

end

end

%% =====================================================================
% ROTATION MATRIX
% =====================================================================

function R = ...
    rotation_matrix_from_a_to_b( ...
    a, ...
    b)

a = ...
    a(:)/norm(a);

b = ...
    b(:)/norm(b);

v = ...
    cross(a,b);

s = ...
    norm(v);

c = ...
    dot(a,b);

if s < 1e-14

    if c > 0

        R = ...
            eye(3);

    else

        if abs(a(1)) < 0.9

            tmp = ...
                [1;0;0];

        else

            tmp = ...
                [0;1;0];

        end

        v = ...
            cross(a,tmp);

        v = ...
            v/norm(v);

        K = ...
            [ ...
             0    -v(3)  v(2); ...
             v(3)  0    -v(1); ...
            -v(2)  v(1)  0];

        R = ...
            eye(3) + ...
            2*K*K;

    end

    return;

end

K = ...
    [ ...
     0    -v(3)  v(2); ...
     v(3)  0    -v(1); ...
    -v(2)  v(1)  0];

R = ...
    eye(3) + ...
    K + ...
    K*K*((1-c)/(s^2));

end

%% =====================================================================
% RANDOM UNIT VECTOR
% =====================================================================

function v = ...
    rand_unit_vector()

x = ...
    randn(1,3);

v = ...
    x/norm(x);

end

%% =====================================================================
% INITIAL RELEASE POINTS
% =====================================================================

function pts = ...
    sample_release_points( ...
    N, ...
    xStart, ...
    releaseRad)

pts = ...
    zeros(N,3);

for i = 1:N

    theta = ...
        2*pi*rand;

    rr = ...
        releaseRad*sqrt(rand);

    y0 = ...
        rr*cos(theta);

    z0 = ...
        rr*sin(theta);

    pts(i,:) = ...
        [xStart y0 z0];

end

end

%% =====================================================================
% FORMAT 2D AXES
% =====================================================================

function format_axes_xy( ...
    xPlotMin, ...
    xPlotMax, ...
    yzPlot)

set(gca, ...
    'FontName','Arial', ...
    'FontSize',14, ...
    'LineWidth',1.0, ...
    'TickDir','out', ...
    'Box','on');

xlabel( ...
    'space, {\itx} (cm)', ...
    'Interpreter','tex', ...
    'FontName','Arial', ...
    'FontSize',16);

ylabel( ...
    'space, {\ity} (cm)', ...
    'Interpreter','tex', ...
    'FontName','Arial', ...
    'FontSize',16);

axis equal;

xlim( ...
    [xPlotMin xPlotMax]);

ylim( ...
    [-yzPlot yzPlot]);

end

%% =====================================================================
% DRAW PANEL: x-y PROJECTION
% =====================================================================

function draw_panel_xy(W)

plot( ...
    [0 0], ...
    [-W/2 W/2], ...
    'r-', ...
    'LineWidth',3);

end

%% =====================================================================
% DRAW RELEASE REGION: x-y PROJECTION
% =====================================================================

function draw_release_segment_xy( ...
    xStart, ...
    releaseRad)

plot( ...
    [xStart xStart], ...
    [-releaseRad releaseRad], ...
    '-', ...
    'Color',[0.45 0.75 0.45], ...
    'LineWidth',3);

end