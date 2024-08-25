function [dates, discounts] = bootstrap(datesSet, ratesSet)
    % bootstrap: bootstraps the discount curve

    % Inputs: 
    % datesSet:     struct containing dates: 
    %               datesSet.depos: vector with deposit expiries
    %               datesSet.futures: 2-col matrix, containing future 
    %               settle dates in the first column, future expiry dates 
    %               in the second
    %               datesSet.swaps: vector containing swaps expiries
    % ratesSet:     struct containing rates: 
    %               ratesSet.depos: matrix with bid-ask deposit rates
    %               ratesSet.futures: matrix with bid-ask future rates
    %               ratesSet.swaps: matrix with bid-ask swap rates

    % Define different year-fraction computation methods
    Act360=2;
    Act365=3;
    Eu30_360=6;

    %Initialize number of futures to use
    numfut=7;
    
    % Initialize dates and discounts
    dates=[datesSet.settlement];
    discounts=[1];

    % Import and save the mid rates
    mid_depos = (ratesSet.depos(:,2) + ratesSet.depos(:,1)) / 2;
    mid_futures = (ratesSet.futures(1:numfut,2) + ratesSet.futures(1:numfut,1)) / 2;
    mid_swaps = (ratesSet.swaps(:,2) + ratesSet.swaps(:,1)) / 2;


    % Build the short end of the curve
    start_futures = find(datesSet.depos > datesSet.futures(1,1), 1,'first');
    dates=[dates; datesSet.depos(1:start_futures-1)];
    discounts = [discounts; 1 ./ (1 + yearfrac(datesSet.settlement, datesSet.depos(1:start_futures-1), Act360) .* mid_depos(1:start_futures-1))];

    % Build the middle end of the curve
    [dates, discounts] = interpolatefutures(datesSet.futures(1:numfut, :), dates, discounts, mid_futures);
   
    % Build the long end of the curve
    
    expiry_swaps=datesSet.swaps;
    zRates = zeroRates(dates, discounts) ./ 100;
    first_zrate = interp1(dates, zRates,expiry_swaps(1), 'linear');
    first_B_swaps =  exp(-first_zrate * yearfrac(dates(1), expiry_swaps(1), Act365));
    
    ind=find(dates>expiry_swaps(1), 1, 'first');
    dates=[dates(1:ind-1); expiry_swaps(1); dates(ind:end)];
    discounts=[discounts(1:ind-1); first_B_swaps; discounts(ind:end)];
    
    BVP = first_B_swaps*yearfrac(datesSet.settlement, expiry_swaps(1), Eu30_360);

    
    for ii = 2:length(expiry_swaps)
    
        B_swaps(ii-1) = (1-mid_swaps(ii)*BVP)/(1+mid_swaps(ii)*yearfrac(expiry_swaps(ii-1), expiry_swaps(ii), Eu30_360));
        BVP = BVP + yearfrac(expiry_swaps(ii-1), expiry_swaps(ii), Eu30_360)*B_swaps(ii-1);
   
    end
    
    dates = [dates; expiry_swaps(2:end)];
    discounts = [discounts; B_swaps'];
end


   