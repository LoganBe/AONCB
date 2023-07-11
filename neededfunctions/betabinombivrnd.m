function [esamp,isamp] = betabinombivrnd(Ke,Ki,nsamp,a0)

% betabinombivrnd(xx,N,nsamp,a0): Generates random samples from the bivariate beta binomial distribution.
%                                 pdf = neCke niCki B(alpha + ke + ki, Ke + Ki - ke - ki + beta)/B(alpha, beta)
%                                 where B is a Beta function
% 
%                                 Under alpha = 0. Uses approximation to
%                                 correctly normalize the pdf given by:
%                                 pdf = neCke niCki B(ke + ki,Ke + Ki - ke - ki + beta)/(psi(beta + Ke + Ki) - psi(beta))
%                                 where psi is the digamma function
%
%                                 Warning: We assume for now that:
%                                 alpha_e == alpha_i
%                                 beta_e == beta_i
%
% Error for approximation for the binomial coef increases for large N.
% Computes binomial coef symbolically when K > 1.5e3, which results in a
% slower computation time.
%
% Input: 
%   Ke -- scalar of total number of excitatory presynaptic neurons (i.e. max number of inputs possible at a given time)
%   Ki -- scalar of total number of excitatory presynaptic neurons (i.e. max number of inputs possible at a given time)
%   nsamp -- Number of samples to take (ns)
%   a0 -- Vector with estimates for alpha and beta given by [alpha, beta]
% Output:
%   esamp -- cell of vector of excitatory samples taken from the bivariate beta binom {number cells each}[nsamp,1]
%   isamp -- cell of vector of inhibtory samples taken from the bivariate beta binom {number cells each}[nsamp,1]
%

assert(Ke > 0 && Ki > 0); %Make sure you have at least 1 synapse for both Ke and Ki

if a0(1) == 0; a0(1) = eps; end

f = @(a,ke,ki) exp(gammaln(Ke+1) + gammaln(Ki+1) + gammaln(ke + ki) + gammaln(a(2)+Ke+Ki-ke-ki)...
                                               -gammaln(ke+1)-gammaln(Ke-ke+1)-gammaln(ki+1)-gammaln(Ki-ki+1)-gammaln(a(2)+Ke+Ki))./(psi(a(2)+Ke+Ki)-psi(a(2)));
xe = 0:Ke; xi = 0:Ki;
[X,Y] = meshgrid(xe,xi);
p = f(a0,X,Y); p = p';
p(1,1) = 0; p(isinf(p))=0;

xy = [];
for i = 1:length(xi)
    xy = [xy; xe', repmat(xi(i),length(xe),1)];
end

esamp = cell(length(nsamp),1);
isamp = cell(length(nsamp),1);
for i = 1:length(nsamp)
    row = randsample(size(xy,1),nsamp(i),true,p(:));
    xy_sample = xy(row,:);
    esamp{i} = xy_sample(:,1); isamp{i} = xy_sample(:,2);
end


