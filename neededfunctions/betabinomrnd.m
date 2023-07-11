function [samples] = betabinomrnd(K,nsamp,a0)

% betabinomrnd(xx,N,nsamp,a0): Generates random samples from the beta binomial distribution.
%                              pdf = nCx B(x + alpha, K - x + beta)/B(alpha, beta)
%                              where B is a Beta function
% 
%                              Under alpha = 0. Uses approximation to
%                              correctly normalize the pdf given by:
%                              pdf = nCx B(x,K - x + beta)/(psi(beta + K) - psi(beta))
%                              where psi is the digamma function
%
% Error for approximation for the binomial coef increases for large N.
% Computes binomial coef symbolically when K > 1.5e3, which results in a
% slower computation time.
%
% Input: 
%   K -- scalar of total number of presynaptic neurons (i.e. max number of inputs possible at a given time)
%   nsamp -- Number of samples to take (ns)
%   a0 -- Vector with estimates for alpha and beta given by [alpha, beta]
% Output:
%   samples -- cell of vector of samples taken from the beta binom {number cells each}[nsamp,1]
%

%Define the betabin pdf
f = @(a,x) exp(gammaln(K+1) + gammaln(a(2)+K-x) - gammaln(a(2)+K) - gammaln(K-x+1)-log(x))./(psi(a(2)+K)-psi(a(2))); 

%Collect samples
samples = cell(length(nsamp),1);
for i = 1:length(nsamp)
    samples{i} = randsample(1:K,nsamp(i),true,f(a0,1:K));
end