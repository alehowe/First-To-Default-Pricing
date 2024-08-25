function zRate=zeroRates(dates,discounts)
% zeroRates: computes the vector of zeroRates at each date

% Inputs:
% dates:        vector containing dates of the curve
% discounts:    vector containing discounts of the curve      

% Compute the zeroRates
Act365=3;
zRate=-log(discounts)./yearfrac(dates(1), dates,Act365).*100;
zRate(1)=0;

end