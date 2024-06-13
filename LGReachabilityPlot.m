if ~exist("solstruct", "var")
    rundir  = 'LGReachability/results/';
    runname = 'K80';   % run to load
    load([rundir runname '/solstruct.mat']);
end
opts = solstruct.opts;

saveplt = false;

% plotting covariance
theta = linspace(0, 2*pi, 100);

% plot x-y trajectory
f = figure(1);
f.Position = [300, 200, 900, 800];
plot(solstruct.opt_m(1,:), solstruct.opt_m(2,:), ...
     '+--', 'Color', 'black', 'DisplayName', '$m_x$'); hold on;
xlabel('x-position [DU]');
ylabel('y-position [DU]');

% plot final solved for covariance
xt   = cos(theta);
yt   = sin(theta);
Srr  = sqrt(10.60)*chol(solstruct.opt_P(1:2,1:2,end));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,end)');
plot(epts(:,1), epts(:,2), 'k-', 'LineWidth', 2, 'DisplayName', '$P_{xx,N}$');

% plot final prescribed covariance 
Srr  = sqrt(10.60)*chol(opts.Pf(1:2,1:2));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2, 'DisplayName', '$P_f$');
legend('Location','best'); legend('AutoUpdate','off');

% plot solved for covariance at various times
for k = 1:3:size(solstruct.opt_P,3)
    Srr  = sqrt(10.60)*chol(solstruct.opt_P(1:2,1:2,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,k)');
    
    figure(1);
    plot(epts(:,1), epts(:,2), 'k');
end

% plot no-zones
xl = xlim;
yl = ylim;
[xno, yno] = meshgrid(xl(1):0.01:xl(2), yl(1):0.01:yl(2)); % Create a grid of x and y values
mask = xno + yno >= 12.75 | xno + 0.1*yno >= 8.75; % Define a logical mask for the region
scatter(xno(mask), yno(mask), 5, 'filled', 'MarkerFaceColor', rgb('LightGray'), 'MarkerFaceAlpha', 0.05); % Plot the points satisfying the mask
%xlim([-5, 10]); ylim([-7.5,16]);

if saveplt
    toTikz([rundir runname '/pos_opt_traj.tex']);
end


% plot x-y velocity trajectory
f = figure(2);
f.Position = [300, 200, 900, 800];
plot(solstruct.opt_m(3,:), solstruct.opt_m(4,:), ...
     '+--', 'Color', 'black', 'DisplayName', '$m_x$'); hold on;
xlabel('x-velocity [VU]');
ylabel('y-velocity [VU]');

% plot final solved for covariance
xt   = cos(theta);
yt   = sin(theta);
Srr  = sqrt(10.60)*chol(solstruct.opt_P(3:4,3:4,end));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, solstruct.opt_m(3:4,end)');
plot(epts(:,1), epts(:,2), 'k-', 'LineWidth', 2, 'DisplayName', '$P_{xx,N}$');

% plot final prescribed covariance 
Srr  = sqrt(10.60)*chol(opts.Pf(3:4,3:4));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(3:4)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2, 'DisplayName', '$P_f$');
legend('Location','best'); legend('AutoUpdate','off');

% plot solved for covariance at various times
for k = 1:3:size(solstruct.opt_P,3)
    Srr  = sqrt(10.60)*chol(solstruct.opt_P(3:4,3:4,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(3:4,k)');
    
    figure(2);
    plot(epts(:,1), epts(:,2), 'k');
end


if saveplt
    toTikz([rundir runname '/vel_opt_traj.tex']);
end

