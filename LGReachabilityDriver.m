%clearvars; close all; clc; cvx_clear;

%% parameter setup
rundir  = 'LGReachability/results/';
runname = 'K1_u1-6';
oldname = 'K0_u1-5';

mkdir([rundir runname '/']);

toPlot = true;

% dimensions and nodes
opts.n = 4;
opts.m = 2;
opts.N = 20;
nU     = opts.N*opts.m;
nX     = (opts.N+1)*opts.n;

% time step
opts.dt = 5/opts.N;

% initial mean, final mean, and final covariance
opts.x0 = [2; 4; 3; 2];
opts.xf = [8; 2; 0; 0];
opts.Pf = diag([0.06; 0.06; 0.006; 0.006]);

% LTI dynamics matrix
Fx      = eye(opts.n);
Fx(1,3) = opts.dt;
Fx(2,4) = opts.dt;
opts.Fx = Fx;

% constant control mapping matrix 
Fu      = [0, 0; ...
           0, 0; ...
           opts.dt, 0; ...
           0, opts.dt];
opts.Fu = Fu;

% process noise mapping matrix
Fw      = eye(opts.n);
opts.Fw = Fw;

% constant process noise covariance
Pww      = 1e-5*eye(opts.n);
opts.Pww = Pww;

% constant measurement mapping matrix
H = [1 0 0 0; ...
     0 1 0 0];
opts.H = H;

% constant measurement noise covariance
Pvv = (2/3)^2.*eye(2);
opts.Pvv = Pvv;

% constant linear gain matrix
lambda = sqrt(Pww(1,1))*opts.dt*opts.dt/sqrt(Pvv(1,1));
r      = (4 + lambda - sqrt(8*lambda + lambda*lambda))/4;
alpha  = 1 - r*r;
beta   = 2*(2 - alpha) - 4*sqrt(1 - alpha);
K = 0.5*opts.Pf*H'/(H*0.5*opts.Pf*H' + Pvv);
%K = [alpha*eye(2); beta*eye(2)/opts.dt];
opts.K = repmat(K, 1, 1, opts.N);

% set of linear gains imported (comment out if using constant gain)
%oldload = load([rundir oldname '/solstruct.mat']);
%opts.K  = oldload.solstruct.Kgains;

% covariance evolution matrix
PHI = (eye(opts.n) - K*H)*Fx;

% new process noise cov
R = (eye(opts.n) - K*H)*Pww*(eye(opts.n) - K*H)' + K*Pvv*K';

% block matrix form of dynamics
Fxbar = eye(opts.n);
PHIbar = eye(opts.n);
Fubar = zeros(opts.n, opts.m*opts.N);
Fwbar = zeros(opts.n, opts.n*opts.N);
Rbar  = [];
for k = 1:opts.N
    Fxbar  = [Fxbar; Fx*Fxbar( end-opts.n + 1:end, : )];

    PHI = (eye(opts.n) - opts.K(:,:,k)*H)*Fx;
    R   = (eye(opts.n) - opts.K(:,:,k)*H)*Pww*(eye(opts.n) - opts.K(:,:,k)*H)' + opts.K(:,:,k)*Pvv*opts.K(:,:,k)';
    Rbar = blkdiag(Rbar, R);

    PHIbar = [PHIbar; PHI*PHIbar( end-opts.n + 1:end, : )] ;
    Fubar  = [Fubar; Fx*Fubar( end-opts.n+1:end,1:(k-1)*opts.m ) Fu zeros(opts.n,(opts.N-k)*opts.m)];
    Fwbar  = [Fwbar; Fw*Fwbar( end-opts.n+1:end,1:(k-1)*opts.n ) Fw zeros(opts.n,(opts.N-k)*opts.n)];
end
Pwwbar = kron(Pww, eye(opts.N+1));
Swwbar = chol(Pwwbar)';
%Rbar   = kron(R, eye(opts.N+1));
Rbar   = blkdiag(Rbar, Pww);
Srrbar = chol(Rbar)';

% control norm constraint
opts.umax = 1.6;                        % max CL control magnitude

% final time mapping matrix for state
EN = zeros(opts.n, nX);
EN(:,end-opts.n+1:end) = eye(opts.n);

% path domain constraints
as = [1, 1; ...
      1, 0.1; ...
      0, 0; ...
      0, 0];
bs = [12.75; 8.75];
ps = [(1 - 0.99)*(0.2/5); (1 - 0.99)*(0.1/5)];

%% convex optimization
%cvx_solver_settings('write', 'dump.task.gz');
cvx_begin sdp
    cvx_solver mosek;
    cvx_precision high;
    % declare optimization variables
    variable V(nU);                % open loop controls
    variable Sxx0(opts.n, opts.n); % initial state covariance SRF
    expression U(opts.N);          % holds control magnitude constraint
    expression P(2,opts.N+1);       % holds state domain constraint

    % build control magnitude constraint
    for i = 1:opts.N
        EU = zeros(opts.m, opts.N*opts.m);  
        EU(:, (i-1)*opts.m + 1: i*opts.m) = eye(opts.m);
        U(i) = norm(EU*V);
    end

    % build final state constraints
    mN = EN*(Fxbar*opts.x0 + Fubar*V); % final mean
    C  = opts.Pf - EN*Rbar*EN';
    X  = [C, EN*PHIbar*Sxx0; ...
          Sxx0'*PHIbar'*EN', eye(opts.n)]; % final covariance LMI

    % declare objecive
    maximize( 2*log_det(Sxx0) );

    subject to
        mN == opts.xf;
        X >= 0;
        for j = 1:2
            for k = 1:opts.N+1
                EX = zeros(opts.n, nX);
                EX(:, (k-1)*opts.n + 1: k*opts.n) = eye(opts.n);
                P(j,k) = as(:,j)'*EX*(Fxbar*opts.x0 + Fubar*V) + norm( [Sxx0'*PHIbar'; Srrbar']*EX'*as(:,j) )*norminv(1 - ps(j));
                P(j,k) <= bs(j);
            end
        end
        U <= opts.umax;
    
cvx_end

cvx_infs = false;
if strcmp(cvx_status, 'Infeasible')
    disp('INFEASIBLE --- NOT STORING OR PLOTTING');
    cvx_infs = true;
end

% collect data
Pxx0 = Sxx0*Sxx0';
opts.S0 = Sxx0;
opts.P0 = Pxx0;

% decompose block matrices
u_OL_arr   = reshape(V, opts.m, []);
opts.uNorm = vecnorm(u_OL_arr, 2, 1);

% solution structure for saving and plotting
solstruct.opt_m  = zeros(opts.n, opts.N+1);
solstruct.opt_P  = zeros(opts.n, opts.n, opts.N+1);
solstruct.u_OL   = u_OL_arr;
solstruct.opts   = opts;

mxbar  = Fxbar*opts.x0 + Fubar*V;
Pxxbar = PHIbar*Pxx0*PHIbar' + Rbar;

solstruct.opt_m(:,1)   = opts.x0;
solstruct.opt_P(:,:,1) = Pxx0;
mkm1 = opts.x0;
Pkm1 = opts.P0;
for k = 2:opts.N+1
    % propagate mean and cov
    mkm = Fx*mkm1 + Fu*solstruct.u_OL(:,k-1);
    Pkm = Fx*Pkm1*Fx' + Pww;
    Pkp = (eye(opts.n) - opts.K(:,:,k-1)*H)*Pkm*(eye(opts.n) - opts.K(:,:,k-1)*H)' + opts.K(:,:,k-1)*Pvv*opts.K(:,:,k-1)';

    % save mean and cov
    solstruct.opt_m(:,k)   = mkm;
    solstruct.opt_P(:,:,k) = Pkp;

    % recursion
    mkm1 = mkm;
    Pkm1 = Pkp;
end

% one more loop to get Kalman gain through time
Pkm1 = opts.P0;
solstruct.Kgains = zeros(opts.n, size(opts.Pvv,1), opts.N);
for k = 1:opts.N
    Pkm = Fx*Pkm1*Fx' + Pww;
    K   = Pkm*H'/(H*Pkm*H' + Pvv);
    solstruct.Kgains(:,:,k) = K;
    Pkp = (eye(opts.n) - K*H)*Pkm*(eye(opts.n) - K*H)' + K*Pvv*K';
    Pkm1 = Pkp;
end

if ~cvx_infs
    save( [rundir runname '/solstruct.mat'], "solstruct" );
    if toPlot
        LGReachabilityPlot;
    end
end













