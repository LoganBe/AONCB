function [ysum] = makepoissinput(N,K,r,t,dt)

% makepoissoninput(N,K,r,t,dt): Creates *uncorrelated* Poissonian inputs with
%                               mean rate r from a population of K neurons that
%                               feed into N neurons. Can use homo- or
%                               heterogenous rates
%                     
%
% Input: 
%   N -- Number of neurons to recieve input [1x1]
%   K -- Number of synapses onto each neuron [1x1]
%   r -- Vector of firing rate for each input cell [1x4]
%   t -- Time vector 
%   dt -- Time step 
%
% Output:
%   ysum = total inputs (sumed spike trains)
%
%

ysum = zeros(length(t),N);
for i = 1:N    
    if length(r) == 1 %IF homogenous rate
        ysum(:,i) = poissrnd(r*dt*1e-3*K,length(t),1); %Sample from poiss dist]
        
    else              %IF heterogenous rate
        ysum(:,i) = poissrnd(r(:,i)*dt*1e-3*K,length(t),1);
    end
end

