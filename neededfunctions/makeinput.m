function [yfe, yfi, ye, yi] = makeinput(t, dt, N ,inputinfo)

% makeffinput(t,dt,N,inputinfo) Create external spiking input.
%                               Input is summed poisson (or compound poisson) process
%                               over multiple input cells to same cell type. 
%                               Input is then filtered for synaptic integration
%
% Input: 
%   t -- time vector [1xn] (ms)
%   dt -- time step (ms)
%   N -- scalar number of cells
%   inputinfo -- strcture with input paramters
%           K.ee = number of excitatory synapses
%           r.ee = rate of excitatory neurons (Hz)
%           K.ei = number of inhibitory synapses
%           r.ei = rate of inhibitoru neurons (Hz)
%           filt = spike train filter type: 'box' or 'exp' (default)
%           tausyn = synaptic time constant (ms)
%           syndel = synaptic time delay (bin): default = 0
%           corrinfo.corr = string if input is correlated: "corr" or "uncorr" -> (default)
%           corrinfo.ae1 = excitatory correlation alpha parameter
%           corrinfo.ae2 = excitatory correlation beta parameter
%           corrinfo.ai1 = inhibitory correlation alpha parameter
%           corrinfo.ai2 = inhibitory correlation beta parameter
%           corrinfo.corridx = crosscorrelation [0 = no cross, 1 = rhoe = rhoi = rhoei]
%
% Output:
%   yfe -- summed spike trains for the [Nee; Nie] excitatory inputs filtered of size [N,length(t)]
%   yfi -- summed spike trains for the [Nei; Nii] inhibitory inputs filtered of size [N,length(t)]
%   ye -- Raw summed excitatory Spike Trains
%   yi -- Raw summed inhibitory Spike Trains

%Check Conditions
if ~isfield(inputinfo,'corr'); inputinfo.corr = 'uncorr'; end %Correlation or not (DEFAULT = NO CORR)
if ~isfield(inputinfo,'filt'); inputinfo.filt = 'box'; end %Synaptic filter type (DEFAULT = BOX)
if ~isfield(inputinfo,'syndel'); syndel = 0; else; syndel = inputinfo.syndel; end %Synaptic delay (DEFAULT + NO DELAY)

%Load In Paramters 
tausyn = inputinfo.tausyn;

%Number of cells
Ke = inputinfo.K.ee; % Excitatory Synapses
Ki = inputinfo.K.ei; % Inhibitory Synapses
K = [Ke, Ki];

%Rates for poisson process. If rate = 0, change flag
rate_e = inputinfo.r.ee; rates{1} = rate_e; % Excitatory Rate
rate_i = inputinfo.r.ei; rates{2} = rate_i; % Inhibitory Rate

%Set rate to 0 if no synaptic inputs
for i = 1:length(K)
    if K(i) == 0; rates{i} = zeros(size(rates{i})); end
end

%Set adjustable rates
be = rates{1}; bi = rates{2};

%If no inputs -> set all outputs to 0 and exit (saves time)
if all(rates{1}(:)==0) && all(rates{2}(:) == 0)
    yfe = zeros(length(t),N);
    yfi = zeros(length(t),N);
    ye = zeros(length(t),N);
    yi = zeros(length(t),N);
    return
end

%Make input (Uncorr = Poisson, Corr = Compound Poisson)
switch inputinfo.corrinfo.corr
    case "uncorr"
        ye = makepoissinput(N,Ke,be,t,dt); %E->E
        yi = makepoissinput(N,Ki,bi,t,dt); %I->E
    case "corr"
        [ye,yi] = makecorrinput(N,Ke,Ki,be,bi,t,dt,inputinfo); %Excitatory
end

% Filter Synaptic Inputs 
if strcmp(inputinfo.filt,'exp') %Exp filter
    yfe = zeros(length(t)+1,N); yfi = zeros(length(t)+1,N); 
    for tt = 2:length(t)+1
        if tt > syndel %Only add filter if we are at time > delay
            yfe(tt,:) = yfe(tt-1,:) - yfe(tt-1,:)/tausyn*dt + ye(tt-1-syndel,:);
            yfi(tt,:) = yfi(tt-1,:) - yfi(tt-1,:)/tausyn*dt + yi(tt-1-syndel,:);
        end
    end 
    yfe(1,:) = []; yfi(1,:) = []; 

elseif strcmp(inputinfo.filt,'box') %Box filter
    yfilt = ones(floor(tausyn/dt),1);
    yfe = filter(yfilt,1,ye);
    yfi = filter(yfilt,1,yi);
end




