function [V, y, oi] = cifv_subthresh(t, dt, N, input, theta, dispflag)

% cifv(t,Ne,Ni,input,theta,dispflag): Simulates voltage from the
% standard conductance based integrate and fire neuron given by:
%
% tau*dV/dt = (Vl-V) + ge(Ve-V) + gi(Vi-V) + I              (1)
%
% Spikes can be generated either through a Rate based description given by
% a nonlinear conditional intensity of lam = k * V(V>0)^n or through a
% threshold method with threshold Vth. 
%
% For subthreshold analysis, use rate based model with reccurent weights of
% 0
%
% Input: 
%   t -- time vector with constant sampling frequency [1xn] (ms)
%   Ne -- scalar for number of excitatory neurons in network
%   Ni -- scalar for number of inhibitory neuron in network
%   input -- strcture giving both excitatory feedfoward inputs
%            (feedfowarde) and inhibitory feedfoward inputs (feedfowardi)
%   theta -- strcture with necessary model paratmers:
%            tauleak -- leak time constant (ms)
%            tausyn -- synaptic time constant (ms)
%            trim -- time index for removel to be in asymtotic state
%            ereversal -- excitatory reversal potential
%            ireversal -- inhibitory reversal potential
%            lreversal -- leak reversal potential
%            ve0 -- initial excitatory voltage
%            vi0 -- initial inhibitory voltage
%            coef -- scaling paramter for conditional intensity (for rate)
%            power -- nonlinear parameter for conditional intensity (for rate)
%            vth -- spiking threshold (for threshold model)
%            ref -- refractory period (ms)
%            vreset -- reset voltage for threshold model
%            resettype -- 'set' (constant value) or 'average' (running
%                         average)
%            modeltype -- 'rate' or 'thresh'
%
% Output:
%   V -- voltage over time given by equation 1 [mxn]
%   y -- spike events over time [mxn]
%   oi -- structure with other output information
%            rate = time varying rate
%            netinput = Ie + Ii
%            a = reccurent filteredspike info
%            ff.e = excitatory feedfoward info
%            ff.i = inhibitory feedfoward info
%            r.e = excitatory recurrent info
%            r.i = inhibitory recurrent info
%            g.e = total excitatory conductance (ff.e + r.e)
%            g.i = total inhibitory conductance (ff.i + r.i)
%            I.e = total excitatory current (g.e.*drivingforce)
%            I.i = total inhibitory current (g.i.*drivingforce)
%
%
 
% Set parameters
tauleak = theta.tauleak; %Leak Time Constant (vectorize for different leaks)
if length(tauleak) > 1; tauleak = tauleak(1)*ones(1,N);
else; tauleak = tauleak*ones(1,N);
end

trimidx = theta.trim; %Time trim idx
Ve_r = theta.ereversal; %Excitatory Reversal Potential
Vi_r = theta.ireversal; %Inhibitory Reversal Potential
if ~isfield(theta,'leak'); Vl = 0; else; Vl = theta.lreversal; end %Leak Voltage (default = 0)
if ~isfield(theta,'ve0'); ve0 = zeros(1,N); else; ve0 = theta.ve0; end %Default inital exciatatory v = 0
if ~isfield(input,'I'); h = 0; else; h = input.I; end %External Input
if nargin < 6; dispflag = 1; end %Used to show model information

%Feedfoward Inputs
g_e = input.feedfowarde; %Exciatory FF Input (1:Ne = EE, Ne+1:end = IE)
g_i = input.feedfowardi; %Inhibitory FF Input (1:Ne = EI, Ne+1:end = II)

%PRINT INFORMATION
if dispflag
    fprintf("Simulating Voltage For:\n")
    fprintf(num2str(N) + " Neurons" + "\n")
    fprintf("Feedfoward Rates (Hz): \n r_ee = " + num2str(theta.r.ee) + " | r_ei = " + num2str(theta.r.ei) + "\n")
end

%Allocation
V = zeros(length(t),N); % Voltage matrix
y = false(length(t),N); %Spike train matrix
netinput = zeros(length(t),N); %Network maxtrix
leak = zeros(length(t),N); %leak matrix

%Total currents
I_e = zeros(length(t),N); %Excitatory
I_i = zeros(length(t),N); %Inhibitory

%Set initial values
V(1,:) = ve0;

%Simulate V
for i = 2:length(t)    
    %driving force
    drive_fe = Ve_r-V(i-1,:); 
    drive_fi = Vi_r-V(i-1,:);
    
    %total current
    I_e(i,:) = g_e(i,:).*drive_fe;
    I_i(i,:) = g_i(i,:).*drive_fi;
    
    %total input
    netinput(i,:) = I_e(i,:) + I_i(i,:);
   
    leak(i,:) = Vl-V(i-1,:);
    
    %Voltage
    V(i,:) = V(i-1,:) + (1./tauleak).*(leak(i,:) + netinput(i,:) + h)*dt; %Voltage
end

%If extended (do to refractory). Remove extra
if size(V,1) > length(t)
    V(length(t)+1:end,:) = [];
end

%Remove trim info for asymtotic stats
V(1:trimidx,:) = []; 
y(1:trimidx,:) = [];

%Save extra info
netinput(1:trimidx,:) = []; oi.netinput = netinput;
leak(1:trimidx,:) = []; oi.leak = leak;

g_e(1:trimidx,:) = []; oi.ff.e = g_e;
g_i(1:trimidx,:) = []; oi.ff.i = g_i;

I_e(1:trimidx,:) = []; oi.I.e = I_e;
I_i(1:trimidx,:) = []; oi.I.i = I_i; 







