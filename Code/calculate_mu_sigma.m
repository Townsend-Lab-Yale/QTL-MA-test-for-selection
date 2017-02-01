% ESTIMATING THE EFFECT DISTRIBUTION OF MUTATIONS FROM MUTATION
% ACCUMULATION DATA 
% Daniel P. Rice and Jeffrey P. Townsend
% 21 November 2011
%  
% Impliments a hill climbing algorithm to find the maximum likelihood
% estimates for the mean and standard deviation of the distribution of
% effect sizes of spontaneous mutations (mu and sig) and the mutation rate.
% 
% INPUT:
% changes is a vector containing the changes in the trait between
% 	measurements in each line
% gens is a vector containing the number of generations between each
%   measurement in each line
% u is the per-generation mutation rate of the trait
% mu0 and sig0 are initial guesses for mu and sig
% mustep0 and sigstep0 are initial step sizes for the algorithm
% maxiter is the number of iterations to run the hill-climbing algorithm
% 
% OUTPUT:
% mu is the mean of the distribution of mutation effects.
% sig is the standard deviation of the distribution of mutation effects.
% 
function [mu, sig, u, likelihood] = estimate_mu_sigma(changes, gens, u0,ustep0, mu0, sig0, mustep0, sigstep0, maxiter)
% Initialize mu, sig and the step sizes
mu = mu0;
sig = sig0;
u = u0;
ustep = ustep0;
mustep = mustep0;
sigstep = sigstep0;


% Keeps track of the acceptance rate over the last 100 proposed steps
last100 = 0;

% Set the maximum number of mutations to consider when calculating the
% likelihood, assuming that the probability of more than cap mutations is
% 0.
cap = capcalculator(u*max(gens));

% Set the factor by which to multiply the step sizes when acceptance rate is low
coolingfactor = .1;

% Calculate the likelihood of the initial guesses for mu and sig
likelihood = mutefflikelihood(mu, sig, changes, gens, u, cap);

i = 0;
while i < maxiter;
    %Propose a step in both parameters and calculate the likelihood of the
    %new values
    newmu = mu + 2*mustep*rand()- mustep;
    newu = u + 2*ustep*rand() - ustep;
    newsig = sig + 2*sigstep*rand() - sigstep;
    newcap = capcalculator(u*max(gens));
    newlikelihood = mutefflikelihood(newmu, newsig, changes, gens, newu, cap);
    
    %Compare new and old likelihoods, if new is greater than old, accept step and
    %increase the acceptance count
    if newlikelihood > likelihood
        mu = newmu;
        sig = newsig;
        cap = newcap;
        u = newu;
        likelihood = newlikelihood;
        last100 = last100 + 1;
    end
    
    i = i + 1;
    
    %Every 100 iterations, display the acceptance rate over the last 100
    %iterations.  If this rate is less than .15, multiply the step sizes by
    %the cooling factor.  Reset the acceptance count
    if mod(i, 100) ==0,
        disp('Accpetance in last 100 iters = ')
        disp(last100/100)
        if last100 < 15
            mustep = mustep*coolingfactor;
            sigstep = sigstep*coolingfactor;
        end
        last100 = 0;
    end
end

% Calculates the cap so that the total probability that is missed by
% approximating the infinite sum is less than .001
function cap = capcalculator(lambda)
cap = 5;
atzero = exp(-lambda);
while (poisscdf(cap,lambda)-atzero)/(1-atzero)<.999
    cap = cap + 1;
end

% Calculates the likelihood of mu and sigma given the MA data (see text)
function l = mutefflikelihood(mu, sig, changes, gens, u, cap)
l = 1;
for i = 1:length(changes)
    l = l*muteffdist(changes(i), u*gens(i), mu, sig, cap);
end

function val = muteffdist(x, lam, mu, sig, cap)
Vals = zeros(1,cap);
for k = 1:cap;
    Vals(k) = poisspdf(k, lam)*normpdf(x, mu*k, sig*sqrt(k));
end
val = sum(Vals);