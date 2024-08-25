function [tau, index]=f_tau(dates_aux, survProbs, intensities,u)
% f_tau: given u, returns the time of default by reversing the probability curve  

% Inputs:
% dates_aux:            dates with the CDSspread payments with initial date
%                       in first position
% survProbs:            vector containing the survival probabilities
% intensities:          vector containing the hazard rates 
% u:                    simulated probability

% Outputs:
% tau:                  time of default
% index:                first index s.t. the survival probabilities are
%                       greater or equal to the simulated u. Also
%                       probablity 1 for the settlement date is considered.
                      

% Define the parameters
Act365=3;

% Define an auxilary vector
survProbs=[1; survProbs];

% Find the probabilities and intensities before u
indexvect=find(survProbs>=u);

% Define auxilary vectors
intensities=[0; intensities];
yearfracVector=yearfrac(dates_aux(1:end-1),dates_aux(2:end),Act365);
yearfracVector_aux=[0;yearfracVector]; 

if (indexvect(end)>=length(dates_aux))
    % In this case there is no deafult, so set the time of default after
    % maturity
    tau=indexvect(end);
else

    dt=yearfrac(dates_aux(1),dates_aux(indexvect(end)),Act365);

    f=@(x) (exp((sum(-intensities(1:indexvect(end)).*yearfracVector_aux(1:indexvect(end)))...
        -intensities(indexvect(end)+1)*(x-dt))))-u; 
    tau=fzero(f,0);

end

% Return the first index s.t. the survival probability is greater or equal
% than u
index=ceil(tau);

end

