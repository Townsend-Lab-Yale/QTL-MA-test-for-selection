% ESTIMATING THE EFFECT DISTRIBUTION OF MUTATIONS FROM MUTATION
% ACCUMULATION DATA 
% Briton Park and Jeffrey P. Townsend
%  
% Implements a hill climbing algorithm to find the maximum likelihood
% estimates for the submodal and supermodal exponential decay parameters of 
% the distribution of effect sizes of spontaneous mutations and the 
% mutation rate (u). 
% 
% 
% INPUT:
% changes is a vector containing the changes in the trait between
% 	measurements in each line
% gens is a vector containing the number of generations between each
%   measurement in each line
% u0 is the initial guess for u
% thetaL0 and thetaR0 are initial guesses for thetaL and thetaR
% thetaLstep0, thetaRstep0, and ustep0 are initial step sizes for the algorithm
% maxiter is the number of iterations to run the hill-climbing algorithm
% 
% OUTPUT:
% thetaL is submodal exponential decay parameter
% thetaR is supermodal exponential decay parameter
% u is the per-generation mutation rate of the trait
% 
function [thetaL, thetaR,u, likelihood] = calculate_LR(changes, gens, u0, ustep0, thetaL0,thetaR0, thetaLstep0,thetaRstep0, maxiter)

% Initialize thetaL, thetaR, u, and the step sizes
thetaL = thetaL0;
thetaR = thetaR0;
u = u0;

thetaLstep = thetaLstep0;
thetaRstep = thetaRstep0;
ustep = ustep0;


% Keeps track of the acceptance rates over the last 100 proposed steps
last100L = 0;
last100R = 0;
last100u = 0;

% Set the maximum number of mutations to consider when calculating the
% likelihood, assuming that the probability of more than cap mutations is
% 0
cap = capcalculator(u*max(gens));

% Set the factor by which to multiply the step sizes when acceptance rate
% is low
coolingfactor = .1;

% Calculate the likelihood of the initial guesses
likelihood = mutefflikelihood(thetaL,thetaR, changes, gens,u, cap);

i = 0;
while i < maxiter;
  %Propose a step in both parameters and calculate the likelihood of the
  %new values
  newthetaL = thetaL + 2*thetaLstep*rand()- thetaLstep;
  newthetaR = thetaR + 2*thetaRstep*rand()- thetaRstep;
  newu = u + 2*ustep*rand() - ustep;
  newcap = capcalculator(u*max(gens));
  newlikelihood = mutefflikelihood(newthetaL,newthetaR, changes, gens, newu, cap);
    
    %Compare new and old likelihoods, if new is greater than old, accept 
    %step and increase the acceptance count
    if newlikelihood > likelihood
        thetaL = newthetaL;
        thetaR = newthetaR;
        u = newu;
        cap = newcap;
        likelihood = newlikelihood;
        last100L = last100L + 1;
        last100R = last100R + 1;
        last100u = last100u +1;
    end
    
    i = i + 1;
    
    %Every 100 iterations, display the acceptance rates over the last 100
    %iterations.  If this rate is less than .15, multiply the step sizes by
    %the cooling factor.  Reset the acceptance count
    if mod(i, 100) ==0,
        disp('Acceptance in last 100 iters = ')
        disp(last100L/100)
        if last100L < 15
            thetaLstep = thetaLstep*coolingfactor;
        end
        last100L = 0;
        
        disp('Acceptance in last 100 iters = ')
        disp(last100R/100)
        if last100R < 15
            thetaRstep = thetaRstep*coolingfactor;
        end
        last100R = 0;
        
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

% Calculates the likelihood of thetaL, thetaR, and u given the MA data
function [l1, l2 ] = mutefflikelihood(thetaL,thetaR, changes, gens, u, cap)
l1 = 1;
l2 = 1; 
for i = 1:length(changes)
    m = muteffdist(changes(i), u*gens(i), thetaL,thetaR, cap);
    l1 = l1*m(1);
    l2 = l2*m(1);
end

end

function [val] = muteffdist(x, lam, thetaL,thetaR, cap)
Vals = zeros(1,cap);
for k = 1:cap
    y = sumLapl(k, thetaL, thetaR);
    Vals(k) = poisspdf(k, lam)*ksdensity(y,x);
end
val = sum(Vals);
end

% Returns a n-dimensional vector y of summed mutational effects
function [y] = sumLapl( k, thetaL, thetaR)
n = 25000;
n_half = n*0.5;
y = zeros(1,n);

for j = 1:k
    temp = zeros(1,n);
    temp(1:n_half) = -1*exprnd(thetaL, 1, n_half);
    temp((n_half+1):n) = exprnd(thetaR,1,n_half);
    y = y + temp(randperm(n));
end
end
