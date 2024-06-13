clearvars; close all; clc;

rundir  = 'LGReachability/results/';
%runname = 'K1';


% save cost for Nruns
Nruns = 51:300;
cost  = zeros(length(Nruns),1);
for nn = Nruns
    oldname = ['K' num2str(nn-1)];
    runname = ['K' num2str(nn)];

    % run optimization
    mkdir([rundir runname '/']);
    LGReachabilityDriver;

    % save cost if problem was solved (not infeasible)
    if cvx_infs
        continue
    else
        cost(nn) = cvx_optval;
    end

    clearvars -except rundir Nruns cost
end