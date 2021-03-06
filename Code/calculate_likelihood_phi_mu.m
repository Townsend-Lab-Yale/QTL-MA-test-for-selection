% QTL Likelihood Calculator
% Briton Park and Jeffrey P. Townsend
% The source code is released under the GPLv3 license

% This code calculates the likelihood of strength of selection c given QTL
% data and the rate parameters of the distribution of mutation
% effects.
% 
% INPUT:
% C is a vector of the values of c to be examined.
% Q is a vector containing the effect sizes the QTL.
% E is a vector containing the standard errors of the QTL.
% theta is a scalar, the estimated rate parameter
% phi is a scalar, the estimated parameter controlling the submodal and supermodal probability mass
% mu is a scalar, the estimated displacement parameter of the distribution
%
% 
% OUTPUT:
% ltable is a matrix whose i,jth entry corresponds to the contribution of
% the ith QTL to the log likelihood of the jth value of c in the vector C.
% sum(ltable) yields the log likelihood of each value of c.
%
% NOTE:  The function can be made parallel by using matlabpool to ask for a
% number of workers and replacing the first for loop with parfor
% 

function ltable = calculate_likelihood_phi_mu(C, Q, E, theta, phi, mu)

% Initialize the output matrix
ltable = zeros(length(Q), length(C));

% Loop over the the values in C
for j = 1: length(C)

    % Initialize the vector of the contribution of each QTL to the
    %likelihood of c
    lvector = zeros(length(Q),1);
    c = C(j);
    
    % Loop over the QTL
    for i = 1:length(Q)
        q = Q(i);
        e = E(i);
        h = 1;
        
        % Draw 1.0*10^5 mutations from the distribution of substitution
        % effects for strength of selection c.
        muts = draw_muts(c, theta, phi, mu);
        
        % Calculate the average likelihood of these mutations.
        L = log(mean(likelihood(muts, q, e)));

        % Find the maximum likelihood of c over all possible numbers of
        % mutations contributing the QTL
        go = 1;
        while go
            % Add a new set of mutation effects to the current set and
            % calculate the new likelihood
            muts2 = muts + draw_muts(c, theta,phi, mu);
            L2 = log(mean(likelihood(muts2, q, e)));

            % If new likelihood is lower than the old likelihood, break the
            % loop and store the old likelihood value.  Otherwise, repeat
            % the loop.
            if L2 < L
                go = 0;
            else
                L = L2;
                muts = muts2;
                h = h+1;
            end
        end

        lvector(i) = L;
    end
    
    % Store the likelihoods in the jth column of the output matrix
    ltable(:,j) = lvector;
end
end

% This function generates 1.0*10^5 samples from the distribution of
% mutation effects given strength of selection c by markov sampling.
function s = draw_muts(c, theta, phi, mu)
delta = .5;
pdf = @(x) probfix2(c,x, theta, phi, mu);
proppdf = @(x,y) unifpdf(y-x,-delta,delta);
proprnd = @(x) x + rand*2*delta - delta;
nsamples = 100000;
x = mhsample(1000*c,nsamples,'pdf',pdf,'proprnd',proprnd,'proppdf', proppdf);
s = x';
end


% This function takes a strength of selection c, and a vector of mutation
% effects X and returns a vector of the probabilities of fixation (See
% text).
function prob = probfix2(c,x, theta, phi, mu)
n =  2.3*10^6;
if (x<mu)
    probability1 = exppdf(abs(x-mu), theta);
else
    probability1 = exppdf(x-mu, theta);
end

if (x <mu)
    probability1 = phi*probability1;
else
    probability1 = (1-phi)*probability1;
end

if c==0
    probability = 1/(2*n)*ones(1,length(x));
else
    probability = (1-exp(-2*c*x))./(1-exp(-4*c*n*x));
end
prob = probability1*probability;

end

% Takes vectors of simulated effects S, observed effects Q, and observed
% standard errors E and returns a vector of the likelihood of each 
%simulated effect.
function L = likelihood(S, Q, E)
L = exp(-((S - Q).^2)./(2*E.^2))./(E * sqrt(2*pi)); 
end
