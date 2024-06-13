clearvars; close all; clc;

saveplt = false;

% load each solstruct and MC stats
scns       = {'SL2_Halo_2day_dx', 'SL2_Halo_2day_dx_obs'};
Nscenario  = length(scns) ;
solstructs = cell(Nscenario,1);
%MCstats    = cell(Nscenario,1);
rundir     = 'NLReachability/results/';
for i = 1:Nscenario
    runname    = scns{i};
    solstructs{i} = load([rundir runname '/solstruct.mat']);
    %MCstats{i}    = load([rundir runname '/MC/MC_stats.mat']);
end
opts = solstructs{1}.solstruct.opts;

% for covariance plots
theta = linspace(0, 2*pi, 100);
xt    = cos(theta); yt = sin(theta);

% plot initial covariances, final covariance, and no zones
f1          = figure(1);
f1.Position = [700, 200, 900, 800];
% plot initial/final covariances and mean trajectories
Corder = colororder("reef");
for k = 1:Nscenario
    % initial covariances
    Srr = sqrt(10.60)*solstructs{k}.solstruct.opts.S0(1:2,1:2);
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstructs{k}.solstruct.opts.x0(1:2)');

    plot(epts(:,1), epts(:,2), 'Color', Corder(k,:)); hold on;
    
    % mean trajectories
    plot(solstructs{k}.solstruct.opt_m(1,:), solstructs{k}.solstruct.opt_m(2,:), ...
     '+--', 'Color', Corder(k,:));

    % final covariances
    Srr = sqrt(10.60)*chol(solstructs{k}.solstruct.opt_P(1:2,1:2, end));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, opts.xf(1:2)');

    plot(epts(:,1), epts(:,2), 'Color', Corder(k,:), 'LineStyle','--'); hold on;
end


% plot final desired covariance
Srr  = sqrt(10.60)*chol(opts.Pf(1:2,1:2));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'k--', 'LineWidth', 2, 'DisplayName', '$P_f$');
%ylim([0, 8]); xlim([0, 10]);

% plot no-zones
% [xno, yno] = meshgrid(-4:0.01:10, -6:0.01:14); % Create a grid of x and y values
% mask = xno + yno >= 12.75 | xno + 0.1*yno >= 8.75; % Define a logical mask for the region
% scatter(xno(mask), yno(mask), 5, 'filled', 'MarkerFaceColor', rgb('LightGray'), 'MarkerFaceAlpha', 0.1); % Plot the points satisfying the mask
% 
% if saveplt
%     cleanfigure;
%     toTikz([rundir '/pos_reachable_pvv.tex']);
% end
% 
% velocity space plots
% plot initial covariances, final covariance, and no zones
f2          = figure(2);
f2.Position = [900, 200, 900, 800];
% plot initial/final covariances and mean trajectories
Corder = colororder("reef");
for k = 1:Nscenario
    % initial covariances
    Srr = sqrt(10.60)*solstructs{k}.solstruct.opts.S0(4:5,4:5);
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstructs{k}.solstruct.opts.x0(4:5)');

    plot(epts(:,1), epts(:,2), 'Color', Corder(k,:)); hold on;

    % mean trajectories
    plot(solstructs{k}.solstruct.opt_m(4,:), solstructs{k}.solstruct.opt_m(5,:), ...
     '+--', 'Color', Corder(k,:));

    % final covariances
    Srr = sqrt(10.60)*chol(solstructs{k}.solstruct.opt_P(4:5,4:5, end));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, opts.xf(4:5)');

    plot(epts(:,1), epts(:,2), 'Color', Corder(k,:), 'LineStyle','--'); hold on;
end


% plot final desired covariance
Srr  = sqrt(10.60)*chol(opts.Pf(4:5,4:5));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(4:5)');
plot(epts(:,1), epts(:,2), 'k--', 'LineWidth', 2, 'DisplayName', '$P_f$');
% 
% 
% if saveplt
%     cleanfigure;
%     toTikz([rundir '/vel_reachable_pvv.tex']);
% end
% 
% 
% % plot control magnitudes
% f2          = figure(3);
% f2.Position = [1100, 200, 900, 800];
% Corder = colororder("reef");
% for k = 1:Nscenario
%     % plot control norms
%     plot(solstructs{k}.solstruct.opts.uNorm, ...
%      'LineStyle', '-', 'Color', Corder(k,:)); hold on;
% 
%     % plot control bound
%     umax = solstructs{k}.solstruct.opts.umax;
%     yline(umax, 'Color', Corder(k,:), 'LineStyle', '--', 'Alpha', 0.4);
% end
% xlim([1, opts.N]);
% 
% 
% if saveplt
%     cleanfigure;
%     toTikz([rundir '/cntrl_reachable_pvv.tex']);
% end



