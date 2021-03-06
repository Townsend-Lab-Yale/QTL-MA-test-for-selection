% ESTIMATING THE EFFECT DISTRIBUTION OF MUTATIONS FROM MUTATION
% ACCUMULATION DATA 
% Briton Park and Jeffrey P. Townsend
% The source code is released under the GPLv3 license

% Implements a hill climbing algorithm to find the maximum likelihood
% estimates for the exponential decay parameter and displacement parameter
% of the distribution of effect sizes of spontaneous mutations, 
%
% INPUT:
% changes is a vector containing the changes in the trait between
% 	measurements in each line
% gens is a vector containing the number of generations between each
%   measurement in each line
% u0 is the intitial gues for u
% theta0 and mu0 are initial guesses for theta and mu
% thetastep0, mustep0, ustep0 is the initial step sizes for the algorithm
% maxiter is the number of iterations to run the hill-climbing algorithm
% 
% OUTPUT:
% theta is the exponential decay parameter
% u is the per-generation mutation rate of the trait% 
% 
function [theta, u, mu, likelihood] = calculate_theta_mu(changes, gens, u0,ustep0, theta0, mu0, thetastep0,mustep0, maxiter)

% Initialize theta, mu, and u and the step sizes
theta = theta0;
mu = mu0;
u = u0;
thetastep = thetastep0;
mustep = mustep0;
ustep = ustep0;

% Keeps track of the acceptance rates over the last 100 proposed steps
last100 = 0;
last100M = 0;
last100u = 0;

% Set the maximum number of mutations to consider when calculating the
% likelihood, assuming that the probability of more than cap mutations is
% 0.
cap = capcalculator(u*max(gens));

% Set the factor by which to multiply the step sizes when acceptance rate
% is low
coolingfactor = .1;

% Calculate the likelihood of the initial guesses for theta, and mu
likelihood = mutefflikelihood(theta, mu, changes, gens, u, cap);

i = 0;
while i < maxiter;
  %Propose a step in all parameters and calculate the likelihood of the
  %new values
  newtheta = theta + 2*thetastep*rand()- thetastep;
  newmu = mu + 2*mustep*rand()- mustep;
  newu = u + 2*ustep*rand() - ustep;
  newcap = capcalculator(newu*max(gens));
  newlikelihood = mutefflikelihood(newtheta, newmu, changes, gens, newu, cap);
    
    %Compare new and old likelihoods, if new is greater than old, accept 
    %step and increase the acceptance count
    if newlikelihood > likelihood
        theta = newtheta;
        mu = newmu;
        u = newu;
        cap = newcap;
        likelihood = newlikelihood;
        last100 = last100 + 1;
        last100M = last100M + 1;
        last100u = last100u + 1;
    end
    
    i = i + 1;
    
    %Every 100 iterations, display the acceptance rates over the last 100
    %iterations.  If this rate is less than .15, multiply the step sizes by
    %the cooling factor.  Reset the acceptance count
    if mod(i, 100) ==0,
        disp('Acceptance in last 100 iters = ')
        disp(last100/100)
        if last100 < 15
            thetastep = thetastep*coolingfactor;
        end
        last100 = 0;
        
        disp('Acceptance in last 100 iters = ')
        disp(last100M/100)
        if last100M < 15
            mustep = mustep*coolingfactor;
        end
        last100M = 0;
        
        disp('Acceptance in last 100 iters = ')
        disp(last100u/100)
        if last100u < 15
            ustep = ustep*coolingfactor;
        end
        last100u = 0;
    end
  
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
end

% Calculates the likelihood of theta and mu given the MA data
function [l1, l2 ] = mutefflikelihood(theta, mu, changes, gens, u, cap)
l1 = 1;
l2 = 1; 
for i = 1:length(changes)
    m = muteffdist(changes(i), u*gens(i), theta, mu, cap);
    l1 = l1*m(1);
    l2 = l2*m(1);
end

end

function [val] = muteffdist(x, lam, theta,mu, cap)
val = zeros(1,cap);
for k = 1:cap
    val(k) = poisspdf(k, lam)*sumLapl(x-k*mu, k, theta);
end
val = sum(val);
end

% Returns the probability of observed mutational effects 
function [prob] = sumLapl(x, n, theta)
prob = 0;
for j = 0:(n-1)
    prob = prob + abs(x)^(n-j-1)*theta^(n-j)*factorial(n+j-1)/(factorial(j)*factorial(n-1)*2^(n+j)*factorial(n-j-1));
end
prob = exp(-1*theta*abs(x))*prob;
end
