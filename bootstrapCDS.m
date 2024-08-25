function [datesCDS, survProbs, intensities] =bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery) 
% bootstrapCDS: bootstraps the survival probabilities and intensities curves
% from a complete set of CDS 

% Inputs:
% datesDF:              vector containing dates of the curve
% discounts:            vector containing discounts of the curve      
% datesCDS              vector containing dates of CDS payments
% spreadsCDS:           vector containing the spreads of the CDS
% flag:                 1 neglecting the accrual,
%                       2 considering the accrual,
%                       3 Jarrow-Turnbull approximation
% recovery:             recovery 

% Outputs:
% datesCDS:             vector containing dates of CDS payments
% survProbs:            vector containing the survival probabilities     
% intensities:          vector containing the intensities

% Parameters
Eu_30_360=6;
Act365=3;

% Find the discounts in the CDS payment dates
discounts=find_discount(datesDF, discounts, datesCDS);

% Define a complete auxilary dates set
datesCDS_complete=[datesDF(1);datesCDS];

% Compute the year fractions between consecutive dates
yearfracVector=yearfrac(datesCDS_complete(1:end-1),datesCDS_complete(2:end),Eu_30_360);

% Computation of the survival probabilities
switch flag
    case {1,2}
        % Build the probabilities vector with an iterative method
        PPrevious=[1];
        for i=1:length(datesCDS)
            [P] = Psolve(discounts(1:i),recovery,spreadsCDS(i),yearfracVector(1:i),PPrevious,flag);
            PPrevious=[PPrevious;P]; 
        end

        survProbs=PPrevious; 

        % Compute the intensities
        intensities=-(1./yearfrac(datesCDS_complete(1:end-1),datesCDS_complete(2:end),Act365)).*log(survProbs(2:end)./survProbs(1:end-1));
    
        % Discard the first probability (1) from the probability curve
        survProbs=survProbs(2:end); 

    case 3
        % Compute the intensities
        intensities=spreadsCDS(1:end)./(1-recovery); 

        % Initialize the probabilities
        survProbs=zeros(length(intensities),1);
        
        % Compute the survival probabilities from the intensities
        for i=1:length(intensities)
            survProbs(i)=exp(-intensities(1:i)'*yearfrac(datesCDS_complete(1:i),datesCDS_complete(2:i+1),Act365)); 
        end
end

end


