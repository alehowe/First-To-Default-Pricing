function [spreadmean, confidence_interval, elapsed_time]=simulateFTD(dates, discounts, datesCDS, CDS_spreads, rho, R, indexmin, M, seed, alpha)
% simulateFTD: simulates a First To Default, finding the spreads s.t. the 
% NPV is zero

% Inputs:
% dates:                vector containing dates of the curve
% discounts:            vector containing discounts of the curve      
% datesCDS              vector containing dates of CDS payments
% CDS_spreads:          vector containing the spreads of the CDS
% rho:                  correlation 
% R:                    recovery vector of each obligor
% expiry:               expiry of the FTD contract
% M:                    number of MC simulations
% seed:                 seed to set for MC simulations
% alpha:                level used for confidence interval (1-alpha)

% Outputs:
% spreadmean:           spread mean of all MC simulations
% confidence_interval:  CI of level alpha for the spread mean
% elapsed_time:         elapsed time of the MC simulations

% Parameters
Act365=3;
Act360=2;
Eu_30_360=6;

% Initialize the number of obligors
num_obligors=length(R);

% Define auxilary dates
dates_aux=[dates(1); datesCDS];
yearfracVector=yearfrac(dates_aux(1:end-1),dates_aux(2:end),Eu_30_360);

% Compute the survival probabilities and intensities. It is done without
% considering the accrual given the very small differences
survProbs=zeros(length(datesCDS), num_obligors);
intensities=zeros(length(datesCDS), num_obligors);
for ii=1:num_obligors
    [datesCDS, survProbs(:,ii), intensities(:,ii)] =bootstrapCDS(dates, discounts, datesCDS, CDS_spreads(:, ii), 2, R(ii));
end

% Initialize the two legs
premium_leg=zeros(1,M);
recovery_leg=zeros(1,M);

% Compute the covariance matrix sigma
sigma=eye(num_obligors)*(1-rho)+rho;

% Compute A s.t. AA'=sigma
A=chol(sigma, 'lower');

% Set the seed
rng(seed);

% Start timing
tic
    
% MC simulations
for jj=1:M
        
        % Sample num_obligors std normal r.v.
        y=randn(1,num_obligors)';
   
        % Compute x and u
        x=A*y;
        u=normcdf(x);
        
        % Compute the intensity and the survival probabilities for each obligor
        tau=zeros(num_obligors, 1);
        index=zeros(num_obligors, 1);
        
        % Simulate taus 
        for ii=1:num_obligors
            [tau(ii), index(ii)]=f_tau(dates_aux, survProbs(:, ii), intensities(:, ii), u(ii));
        end

        % Take the minimum tau, index and the obligor who first defaults
        tau_min=min(tau);
        index_min=min(index);
        first_defaulter=find(tau_min==tau);
      
        if floor(tau_min)>=indexmin
            % If the default is after the maturity of the contract the
            % recovery leg is 0 and the premium leg if the full BPV

            recovery_leg(jj)=0;
            premium_leg(jj)=sum(find_discount(dates, discounts, dates_aux(2:end)).*yearfracVector);
                                   
        else
            % If the default is before maturity the premium leg is up to
            % tau and the recovery leg is calculated in tau

            discount_at_tau=find_discount(dates, discounts, tau_min*365+dates_aux(1));
            recovery_leg(jj) = (1-R(first_defaulter))*discount_at_tau;
            
            if (tau_min<1)
                yearfracts=0;
            else  
                yearfracts=yearfrac(dates_aux(1:index_min-1), dates_aux(2:index_min), Act365);
            end
            premium_leg(jj) = sum(find_discount(dates, discounts, dates_aux(2:index_min)).*yearfracts)+...
                                discount_at_tau*(tau_min-(index_min-1));     
            end
end

% End timing
elapsed_time=toc;

% Mean and for the mean at level alpha
mean_recovery= mean(recovery_leg);
mean_premium=mean(premium_leg); 
spreadmean=mean_recovery/mean_premium;

% Confidence interval computation
cov_ratio=cov(recovery_leg, premium_leg);
dev_std=sqrt(cov_ratio(1,1)/mean_recovery^2+cov_ratio(2,2)/mean_premium^2-2*cov_ratio(1,2)/(mean_premium*mean_recovery));
error=dev_std/sqrt(M);
t_stat=tinv(1-alpha/2, M-1);
confidence_interval=[spreadmean-spreadmean*error*t_stat, spreadmean+spreadmean*error*t_stat];
   
end

