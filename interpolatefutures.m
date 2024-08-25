function [dates, discounts] = interpolatefutures(future_dates, dates, discounts, mid_futures)
    % interpolatefutures: builds the middle end of the curve

    % Inputs: 
    % future_dates: 2-col matric containing future 
    %               settle dates in the first column, future expiry dates 
    %               in the second
    % dates:        vector containing dates of the short end of the
    %               curve (from deposit dates)
    % discounts:    vector containing discounts of the short end of the
    %               curve (from deposit rates)             
    % mid_futures:  mean of bid-ask of futures

    % Define different year-fraction computation methods:
    Act360 = 2;
    Act365 = 3;
    Eu30_360=6;

    % Extract and save dates
    settle_futures = future_dates(:, 1);
    expiry_futures = future_dates(:, 2);
    
    % Calculate zero rates for deposits
    zRates = zeroRates(dates, discounts) ./ 100;

    % Interpolate zero rates for the first future contract
    settle_futures = future_dates(:, 1);
    expiry_futures = future_dates(:, 2);
    

    % Compute discount factors 
    B_end = 1 ./ (1 + yearfrac(settle_futures, expiry_futures, Act360) .* mid_futures);
    
    %Compute discount factor related to the first future
    first_zRate = interp1(dates, zRates, settle_futures(1), 'linear');
    first_B = exp(-first_zRate * yearfrac(dates(1), future_dates(1, 1), Act365));
    
    % Update discounts and dates
    discounts = [discounts; first_B * B_end(1)];
    dates = [dates; expiry_futures(1)];
    zRates = zeroRates(dates, discounts) ./ 100;

    % Interpolate/extrapolate zero rates and compute discount factors for subsequent future contracts
    for ii = 2:size(future_dates, 1)
        ind = settle_futures(ii) > expiry_futures(ii - 1);
        
        switch ind
            case 0
                % linear interpolation
                nth_zRate = interp1(dates, zRates, future_dates(ii, 1), 'linear');
            case 1
                % flat extrapolation 
                nth_zRate = interp1(dates, zRates, future_dates(ii, 1), 'previous', 'extrap');
        end

        nth_B = exp(-nth_zRate * yearfrac(dates(1), future_dates(ii, 1), Act365));

        % Update discounts and dates
        discounts = [discounts; nth_B * B_end(ii)];
        dates = [dates; future_dates(ii, 2)];

        % Update zero rates
        zRates = zeroRates(dates, discounts) ./ 100;
    end

end







