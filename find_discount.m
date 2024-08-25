function discount = find_discount(dates, discounts, searching_dates)
    % discount: finds the discount vector of the corresponding
    % searching_dates

    % Inputs:
    % dates:            vector containing dates of the curve
    % discounts:        vector containing discounts of the curve  
    % searching_dates:  vector of dates corresponding to the discounts
    %                   needed

    % Initialize time set
    Act365 = 3;

    % Preallocate discounts vector
    discount = zeros(size(searching_dates)); 
    
    for i = 1:numel(searching_dates)
        % Find index of searching_date in dates array
        index = find(searching_dates(i) == dates, 1, 'first');

        % If searching_date is found in dates array
        if ~isempty(index)
            discount(i) = discounts(index);
        
        % If searching_date is within the range of dates
        elseif searching_dates(i) < dates(end)
            zero_rates = zeroRates(dates, discounts)/100;
            z_rate = interp1(dates, zero_rates, searching_dates(i), 'linear');
            discount(i) = exp(-z_rate * yearfrac(dates(1), searching_dates(i), Act365));

        % If searching_date is beyond the last date in dates array
        else
            zero_rates = zeroRates(dates, discounts)/100;
            z_rate = interp1(dates, zero_rates, searching_dates(i), 'linear', 'extrap');
            discount(i) = exp(-z_rate * yearfrac(dates(1), searching_dates(i), Act365));
        end
    end
end


