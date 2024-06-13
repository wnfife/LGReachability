clearvars; close all; clc;

% load case solution struct
rundir  = 'LGReachability/results/';
runname = 'K1';
load([rundir runname '/solstruct.mat']);
opts = solstruct.opts;

% make MC directory
mcdir = [rundir runname '/MC/'];
mkdir(mcdir);

rng(100);

% number of samples
Nmc = 20000;

% storage
MC_xOL.all = zeros(opts.N, Nmc, opts.n);
MC_uOL.all = zeros(opts.N-1, Nmc, opts.m);
MC_Pxx     = zeros(opts.n, opts.n, opts.N);

% process noise SRF for sampling noise
Sww = chol(opts.Pww)';
Svv = chol(opts.Pvv)';
% srf for sampling
S0 = (1).*opts.S0;
% start MC loop
for i = 1:Nmc
    % sample from Gaussian
    x0i = opts.x0 + S0*randn(opts.n,1);

    % start initial cov
    Pkm1 = S0*S0';
    MC_Pxx(:,:,1) = Pkm1;

    % start sim loop
    MC_xOL.all(1,i,:) = x0i;
    xkm1 = x0i;
    mkm1 = solstruct.opt_m(:,1);
    for k = 1:opts.N
        % extract control
        u_OL = solstruct.u_OL(:,k);

        % save control
        MC_uOL.all(k,i,:) = u_OL;

        % propagate sample
        wkm1 = Sww*randn(opts.n,1);
        xkm  = opts.Fx*xkm1 + opts.Fu*u_OL;

        % propagate sample filter covariance
        Pkm = opts.Fx*Pkm1*opts.Fx' + opts.Pww;

        % propagate truth
        mk = opts.Fx*mkm1 + opts.Fu*u_OL + wkm1;

        % update sample and covariance
        zk  = opts.H*mk + Svv*randn(2,1);
        K   = Pkm*opts.H'/(opts.H*Pkm*opts.H' + opts.Pvv);
        %K   = opts.K(:,:,k);
        xkp = xkm + K*(zk - opts.H*xkm);
        Pkp = (eye(opts.n) - K*opts.H)*Pkm*(eye(opts.n) - K*opts.H)' + K*opts.Pvv*K';

        % save samples
        MC_xOL.all(k+1, i, :) = xkp;
        MC_Pxx(:,:,k+1) = Pkp;

        % recursion
        xkm1 = xkp;
        mkm1 = mk;
        Pkm1 = Pkp;
    end
end

% save
MC_stats.MC_xOL = MC_xOL;
MC_stats.MC_uOL = MC_uOL;
MC_stats.MC_Pxx = MC_Pxx;
save([mcdir 'MC_stats'], 'MC_stats');

LGReachabilityMCPlot;