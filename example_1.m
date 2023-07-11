% Example of the All-or-none conductance-based (AONCB) neuron under 3 conditions:
% Uncorrelated inputs, Excitatory and Inhibitory correlations (no cross corr)
% and Excitatory + Inhibitory correlatiosn with cross correlations.
% Plots example of the 3 traces along with the theoretical estimates using
% both the full form and the smallweight approximations

% The AONCB neuron is defined by the standard conductance model given by:
%
%               CdV/dt = G(Vl-V) + G_e(Ve-V) + G_i(Vi-V) + I
%
% where the conductances G_e/G_i are given by an all-or-none process where
% given an input, the channel opens up with size given by synaptic weight W 
% for tausyn miliseconds then closes.
% G_e exampel:
%
%        |---| <- tausyn 
%    |   _____
%    | W |   |
%    |___|   |____
%    |____________
%               time
%
% We define small weight estimates of the moments as (see paper for full form):
% mean = (Ke*re*we*Ve + Ki*ri*wi*Vi + I)/(1 + Ke*re*we + Ki*ri*wi)
%
% var = (Ke*re*we^2*(1+rhoe*(Ke-1))*(Ve-mean)^2 + Ki*ri*wi^2*(1+rhoi*(Ki-1))*(Vi-mean)^2 - 2*rhoei*sqrt(re*ri)*Ke*we*Ki*wi*(Ve-mean)*(mean-Vi))/
%       (2*tauleak*(1 + Ke*re*we + Ki*ri*wi))
%
% Paper: Exact analysis of the subthreshold variability for conductance-based neuronal models with synchronous synaptic inputs. https://arxiv.org/abs/2304.09280

rng(2312);
set(0,'DefaultAxesFontSize',20)
%% Set Model Paramters
dt = 0.01; %Time step (ms)
tmax = 2e3; %Max time (ms)
ttrim = 100; %Time to remove (ms) Puts neuron in asymtotic state
t = 0:dt:tmax+ttrim; % time vector
tt = t(1:end-ttrim/dt-1)*1e-3; % Time Vector Trimed in in Seconds (used for plotting)
trimidx = binsearch(t,ttrim); % trim time idx

%Time Constants
tauleak = 15; % Leak Time Constant (ms)
tausyn = 0.02; % Synaptic Input Time Constant (ms)

%Biophysical Paramters (Normalized to Vl = 0)
Vl = 0; % Leak voltage
Ve_r = 60; %Excitatory reversal potential
Vi_r = -10; %Inhibitory reversal potential

%Number of cells (trials)
N = 500;

%Set homogenous synaptic weights (Set so EPSP -> 0.05 mV; IPSP -> 4*EPSP)
Wee = 0.015; % E->E (ms)
Wei = 0.06;  % I->E (ms)
W_f = [Wee,Wei];

%Scaled weight for each synapse
J_f = zeros(N,2);
J_f(:,1) = Wee./tausyn;
J_f(:,2) = Wei./tausyn;

%Excitatory synaptic informaiton
Ke = 1000; %Number Excitatory synapses
re = 10; % (Hz)

%Inhibitory synaptic information
Ki = 1000; %Number Inhibitory Synapses
ri = 10; %(Hz)

%Syn Integration Type (box,exp)
filttype = "box";

% Set Parameters Into Structure
theta.tauleak = tauleak;
theta.tausyn = tausyn;
theta.trim = trimidx;
theta.ereversal = Ve_r;
theta.ireversal = Vi_r;
theta.Wf = W_f;
theta.WeightFeedFoward = J_f;
theta.K.ee = Ke; theta.r.ee = re;
theta.K.ei = Ki; theta.r.ei = ri;
theta.filt = filttype;

%% SIMULATE EXAMPLE WITH NO CORR
corrtype = "uncorr";
theta.corrinfo.corr = corrtype;

%Make feedfoward inputs
[ye_nc, yi_nc] = makeinput(t,dt,N,theta); %[excitatory input spike counts, inhibitory input spike counts]
wye_nc = ye_nc.*J_f(:,1)'; wyi_nc = yi_nc.*J_f(:,2)'; %Weighted inputs

input.feedfowarde = wye_nc;
input.feedfowardi = wyi_nc;
input.I = 0; 

%sim V
Vnc = cifv_subthresh(t, dt, N, input, theta,0);
meanVnc = mean(Vnc(:)); %mean V
varVnc = mean(var(Vnc,[],2)); %var V

[meanVnct, varVnct] = statest_full(theta); %Use full theoretical moment estimate 
[meanVncsw,varVncsw] = statest_sw(theta);  %Use Small weight approximation moment estimate

%% SIMULATE EXAMPLE WITH CORR BUT NO CROSS CORR
corrtype = "corr";
rhoe = 0.03; rhoi = 0.03; %Excitatory and inhibitory marginal correlations
xcorridx = 0; %No Cross Corr

%Sets paramters
theta.corrinfo.corr = corrtype;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-theta.corrinfo.ae1-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rhoi)-theta.corrinfo.ai1-1;
theta.corrinfo.corridx = xcorridx;

%Make feedfoward inputs
[ye_wc, yi_wc] = makeinput(t,dt,N,theta); %[excitatory input spike counts, inhibitory input spike counts]
wye_wc = ye_wc.*J_f(:,1)'; wyi_wc = yi_wc.*J_f(:,2)'; %Weighted inputs

input.feedfowarde = wye_wc;
input.feedfowardi = wyi_wc;
input.I = 0;

%sim V
Vwc = cifv_subthresh(t, dt, N, input, theta);
meanVwc = mean(Vwc(:)); %Mean V
varVwc = mean(var(Vwc,[],2)); %varV

[meanVwct, varVwct] = statest_full(theta); %Use full theoretical moment estimate 
[meanVwcsw,varVwcsw] = statest_sw(theta,rhoe,rhoi,0);  %Use Small weight approximation moment estimate

%% SIMULATE EXAMPLE WITH CORR AND MAX CROSS CORR (NOTE THIS ONLY CURRENTLY WORKS FOR RHOEI = RHOE = RHOI; KE = KI; RE + RI)
corrtype = "corr";
rhoe = 0.03; rhoi = 0.03;%Excitatory and inhibitory marginal correlations
xcorridx = 1; %Set Max Corr -> Implies that rhoei = rhoe = rhoi
rhoei = 0.03;

%Sets paramters
theta.corrinfo.corr = corrtype;
theta.corrinfo.ae1 = 0; theta.corrinfo.ae2 = (1/rhoe)-theta.corrinfo.ae1-1;
theta.corrinfo.ai1 = 0; theta.corrinfo.ai2 = (1/rhoi)-theta.corrinfo.ai1-1;
theta.corrinfo.corridx = xcorridx;

%Make feedfoward inputs
[ye_xc, yi_xc] = makeinput(t,dt,N,theta); %[excitatory input spike counts, inhibitory input spike counts]
wye_xc = ye_xc.*J_f(:,1)'; wyi_xc = yi_xc.*J_f(:,2)'; %Weighted inputs

input.feedfowarde = wye_xc;
input.feedfowardi = wyi_xc;
input.I = 0; %External current injection

%sim V
Vxc = cifv_subthresh(t, dt, N, input, theta);
meanVxc = mean(Vxc(:));
varVxc = mean(var(Vxc,[],2));

[meanVxct, varVxct] = statest_full(theta); %Use full theoretical moment estimate 
[meanVxcsw,varVxcsw] = statest_sw(theta,rhoe,rhoi,rhoei);  %Use Small weight approximation moment estimate

%% Set Table
names = ["No Corr";"Within Corr";"Cross Corr"];
meanVsim = [meanVnc; meanVwc; meanVxc];
meanVt = [meanVnct; meanVwct; meanVxct];
meanVsw = [meanVncsw; meanVwcsw; meanVxcsw];
varVsim = [varVnc; varVwc; varVxc];
varVt = [varVnct; varVwct; varVxct];
varVsw = [varVncsw; varVwcsw; varVxcsw];

T = table(names, meanVsim, meanVt, meanVsw, varVsim, varVt,varVsw);
%% Plot figures
tidx = binsearch(tt,1);
limtemp = [Vnc(:,1); Vwc(:,1); Vxc(:,1)];
figure(1); clf;
subplot(2,8,[1,6]); hold on
plot(tt(1:tidx)*1e3,Vnc(1:tidx,1))
plot(tt(1:tidx)*1e3,Vwc(1:tidx,1))
plot(tt(1:tidx)*1e3,Vxc(1:tidx,1))
xlabel('Time (ms)'); ylabel('Voltage (mV)'); 
legend({'No Corr','Within Corr','Cross Corr'},'box','off')
ylim([min(limtemp(:))+min(limtemp(:))*0.4,max(limtemp(:))+max(limtemp(:))*0.4])

subplot(2,8,[7,8]); hold on
histogram(Vnc(:),'binwidth',0.5,'normalization','pdf','orientation','horizontal')
histogram(Vwc(:),'binwidth',0.5,'normalization','pdf','orientation','horizontal')
histogram(Vxc(:),'binwidth',0.5,'normalization','pdf','orientation','horizontal')
yticks([]); ylim([min(limtemp(:))+min(limtemp(:))*0.4,max(limtemp(:))+max(limtemp(:))*0.4])

subplot(2,8,[9,11.8])
bar([mean(meanVnc),mean(meanVncsw);mean(meanVwc),mean(meanVwcsw);mean(meanVxc),mean(meanVxcsw)]); box off
ylabel('Mean V'); legend({'Sim','Theory'},'box','off'); ylim([0,10]);  xticklabels({'NC','WC','CC'})
subplot(2,8,[13.2,16])
bar([mean(varVnc),mean(varVncsw);mean(varVwc),mean(varVwcsw);mean(varVxc),mean(varVxcsw)]); box off
ylabel('Var V'); xticklabels({'NC','WC','CC'})
