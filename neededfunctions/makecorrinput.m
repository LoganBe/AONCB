function [ye, yi] = makecorrinput(N, Ke, Ki, re, ri, t, dt, theta)

% makecorrinput(N,Ke,Ki,r_e,r_i,t,dt,theta): Generates correlated inputs,
%                                            through a compound poisson process with a beta measure
%                                            with within population correlations (rho_e, rho_i)
%                                            and cross population correlations (rho_ei)
%
% Inputs:
%   N -- Total number of output neurons [1x1]
%   Ke -- Total number of excitatory synapses [1x1]
%   Ki -- Total number of inhibitory synapses [1x1]
%   re -- Excitatory input rate vector [nxm] (Hz)
%   ri -- Inhibitory input rate vector9 [nxm] (Hz)
%   t -- time vector [1xn] (ms)
%   dt -- time step (ms)
%   theta -- structure with parameters given by:
%            dist: Indicates jump size distribution - "betabin" (default),
%            "unif", "exp", "gamma", "poiss", "lognormal","wbl"
%            corrinfo.corridx: level of corr across pop (rho_ei) [0 = no corr, 1 = max corr]
%            corrinfo.ae1: first excitatory paramter (i.e alpha in betabin)
%            corrinfo.ae2: second excitatory paramter (i.e beta in betabin)
%            corrinfo.ai1: first inhibitory paramter (i.e alpha in betabin)
%            corrinfo.ai2: second inhibitory paramter (i.e beta in betabin)
%
% Outputs:
%  ye -- Summed excitatory spike train [mxn]
%  yi -- Summed inhibitory spike train [mxn]
%
%

%Set defaults
if isfield(theta.corrinfo,'corridx'); ci = theta.corrinfo.corridx; else; ci = 0; end %default cross = 0

%If only 1 input type, make sure ci = 0
if all(re(:) == 0) || all(ri(:) == 0); ci = 0; end

%If no inputs at all, quick return zeros
if all(re(:) == 0) && all(ri(:) == 0)
    ye = zeros(length(t),N);
    yi = zeros(length(t),N);
    return
end

%Beta-bin paramters (we set alpha = 0)
alpha_e = theta.corrinfo.ae1; beta_e = theta.corrinfo.ae2;
alpha_i = theta.corrinfo.ai1; beta_i = theta.corrinfo.ai2;

%Adjust rates to match with correlatios
be = re*beta_e*(psi(beta_e+Ke)-psi(beta_e))*1e-3; %adjusted E rate where psi is the digamma function
bi = ri*beta_i*(psi(beta_i+Ki)-psi(beta_i))*1e-3; %adjusted I rate

%Marginal spike time info
spidx_e = be.*(1-ci).*dt > rand(length(t),N); %excitatory E only spike time 
spidx_i = bi.*(1-ci).*dt > rand(length(t),N); %Inhibitory I only spike time

%Joint spike times
if be >= bi
    spidx_ei1 = be.*ci.*(bi./be)*dt > rand(length(t),N); 
    spidx_ei2 = be.*ci.*(1-bi./be)*dt > rand(length(t),N);
else
    spidx_ei1 = bi.*ci.*(be./bi)*dt > rand(length(t),N);
    spidx_ei2 = bi.*ci.*(1-be./bi)*dt > rand(length(t),N);
end

%Total spike trian vector
ye = zeros(length(t),N); 
yi = zeros(length(t),N);

%Marginal Spike Samples (Don't need if using max cross corr)
if ci ~= 1 
    esample = betabinomrnd(Ke,sum(spidx_e),[alpha_e,beta_e]); %sample E only
    isample = betabinomrnd(Ki,sum(spidx_i),[alpha_i,beta_i]); %sample I only
end

%Join spike samples (don't need if no cross corr)
if ci~= 0
    if be >= bi
        eisample = betabinomrnd(Ke,sum(spidx_ei2),[alpha_e,beta_e]); %sample E only marginal
    else 
        eisample = betabinomrnd(Ki,sum(spidx_ei2),[alpha_i,beta_i]); %sample I only margnial
    end 
    [ejsample,ijsample] = betabinombivrnd(Ke,Ki,sum(spidx_ei1),[alpha_e,beta_e]); %sample joint EI
end

%Plug in samples into spike train wwhen we have a spike
for nn = 1:N
    % Indep E and I Spikes 
    if ci ~= 1
        ye(spidx_e(:,nn),nn) = esample{nn}; %Sample P(We) for (1-alpha)%
        yi(spidx_i(:,nn),nn) = isample{nn};  %Sample P(Wi) for (1-alpha)%
    end

    % Joint with indep spikes
    if ci ~= 0
        if be >= bi
            ye(spidx_ei2(:,nn),nn) = eisample{nn};
        else
            yi(spidx_ei2(:,nn),nn) = eisample{nn};
        end
        ye(spidx_ei1(:,nn),nn) = ejsample{nn};
        yi(spidx_ei1(:,nn),nn) = ijsample{nn};
    end
end
        
        
