if ~exist("solstruct", "var")
    rundir  = 'LGReachability/results/';
    runname = 'umax_2_sig2-3';   % run to load
    load([rundir runname '/solstruct.mat']);
    load([rundir runname '/MC/MC_stats.mat']);
end

mcdir = [rundir runname '/MC/'];

opts = solstruct.opts;

saveplt = false;

% plot x-y OL trajectories
f2          = figure(2);
f2.Position = [700, 200, 900, 800];
plot(MC_stats.MC_xOL.all(:,1:100:end,1), MC_stats.MC_xOL.all(:,1:100:end,2), 'Color', [rgb('DarkSlateGray'), 0.2] ); hold on;
xlabel('x-position [DU]'); ylabel('y-position [DU]');

% plotting covariance
theta = linspace(0, 2*pi, 100);
xt    = cos(theta);
yt    = sin(theta);

% plot optimized covariance at various times
for k = 1:2:size(solstruct.opt_P,3)
    Srr  = sqrt(10.60)*chol(solstruct.opt_P(1:2,1:2,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,k)');
    
    figure(2);
    plot(epts(:,1), epts(:,2), 'k');
end

% plot MC covariance at various times
for k = 1:2:size(MC_stats.MC_Pxx,3)
    Srr  = sqrt(10.60)*chol(MC_stats.MC_Pxx(1:2,1:2,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,k)');
    
    figure(2);
    plot(epts(:,1), epts(:,2), 'r');
end

% plot final prescribed covariance 
Srr  = sqrt(10.60)*chol(opts.Pf(1:2,1:2));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2);

% plot final optimized covariance 
Srr  = sqrt(10.60)*chol(solstruct.opt_P(1:2,1:2,end));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'k--', 'LineWidth', 1);

% plot no-zones
xl = xlim;
yl = ylim;
[xno, yno] = meshgrid(xl(1):0.01:xl(2), yl(1):0.01:yl(2)); % Create a grid of x and y values
mask = xno + yno >= 12.75 | xno + 0.1*yno >= 8.75; % Define a logical mask for the region
scatter(xno(mask), yno(mask), 5, 'filled', 'MarkerFaceColor', rgb('LightGray'), 'MarkerFaceAlpha', 0.05); % Plot the points satisfying the mask
%xlim([0, 10]); ylim([0, 10]);

if saveplt
    toTikz([mcdir '/MC_xCL.tex']);
end


% plot x-y velocity OL trajectories
f3          = figure(3);
f3.Position = [900, 200, 900, 800];
plot(MC_stats.MC_xOL.all(:,1:100:end,3), MC_stats.MC_xOL.all(:,1:100:end,4), 'Color', [rgb('DarkSlateGray'), 0.2] ); hold on;
xlabel('x-velocity [VU]'); ylabel('y-velocity [VU]');

% plotting covariance
theta = linspace(0, 2*pi, 100);
xt    = cos(theta);
yt    = sin(theta);

% plot optimized covariance at various times
for k = 1:2:size(solstruct.opt_P,3)
    Srr  = sqrt(10.60)*chol(solstruct.opt_P(3:4,3:4,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(3:4,k)');
    
    figure(3);
    plot(epts(:,1), epts(:,2), 'k');
end

% plot MC covariance at various times
for k = 1:2:size(MC_stats.MC_Pxx,3)
    Srr  = sqrt(10.60)*chol(MC_stats.MC_Pxx(3:4,3:4,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(3:4,k)');
    
    figure(3);
    plot(epts(:,1), epts(:,2), 'r');
end

% plot final prescribed covariance 
Srr  = sqrt(10.60)*chol(opts.Pf(3:4,3:4));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(3:4)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2);

% plot final optimized covariance 
Srr  = sqrt(10.60)*chol(solstruct.opt_P(3:4,3:4,end));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(3:4)');
plot(epts(:,1), epts(:,2), 'k--', 'LineWidth', 1);

if saveplt
    toTikz([mcdir '/MC_VxCL.tex']);
end
