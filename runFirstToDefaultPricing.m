% runproject
% Howe Alessandro John

clear all;
close all;
clc;

%% Settings
formatData='dd/mm/yyyy'; 

%% Read market data
% This fuction works on Windows OS.

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

%% Bootstrap
% dates includes SettlementDate as first date

[dates, discounts]=bootstrap(datesSet, ratesSet); 

%% Complete set of CDS spreads for ISP and UCG

% Spline interpolation for ISP
CDS_fixed_dates_tointerpolate=[datesSet.swaps(1:5); datesSet.swaps(7)];
CDS_spreadsISP_tointerpolate=[29; 32; 35; 39; 40; 41]*10^(-04);
CDS_spreadsISP6y=interp1(CDS_fixed_dates_tointerpolate, CDS_spreadsISP_tointerpolate, datesSet.swaps(6), 'linear');

% Build the CDS spreads for ISP
CDS_fixed_dates=datesSet.swaps(1:7);
CDS_spreads_ISP=[CDS_spreadsISP_tointerpolate(1:5); CDS_spreadsISP6y; CDS_spreadsISP_tointerpolate(end)];

% Spline interpolation for UCG (the dates are the same)
CDS_spreadsUCG_tointerpolate=[34; 39; 45; 46; 47; 47]*10^(-04);
CDS_spreadsUCG6y=interp1(CDS_fixed_dates_tointerpolate, CDS_spreadsUCG_tointerpolate, datesSet.swaps(6), 'linear');

% Build the CDS spreads for UCG
CDS_spreads_UCG=[CDS_spreadsUCG_tointerpolate(1:5); CDS_spreadsUCG6y; CDS_spreadsUCG_tointerpolate(end)];

% Plot the two spread curves
figure;
plot(datetime(CDS_fixed_dates, 'ConvertFrom', 'datenum'), CDS_spreads_ISP*10^4, '-o');
hold on;
plot(datetime(CDS_fixed_dates, 'ConvertFrom', 'datenum'), CDS_spreads_UCG*10^4, '-o');
title("Spread curve (bps)");
legend("Spread ISP", "Spread UCG", location="northwest");
grid on;
hold off;

%% CDS bootstrap

% Parameters
Act365=3;
Act360=2;
Eu_30_360=6;

% Define the two obligors
obligors = ["ISP", "UCG"];

% Define the two recovery values
R_values = [0.4, 0.45];  

% Define the two sets of spreads
CDS_spreads=[CDS_spreads_ISP, CDS_spreads_UCG];


% Compute the intensities and survival probabilities for each obligor
for idx = 1:2
    obligor=obligors(idx);
    R = R_values(idx);
    spreads = CDS_spreads(:, idx);

    % Plot for the current obligor
    figure;

    % Compute and plot neglecting accrual
    flag = 1;
    [datesCDS, survProbs, intensitiesNOaccrual] = bootstrapCDS(dates, discounts, CDS_fixed_dates, spreads, flag, R); 
    fprintf("Intensities neglecting accrual for %s are:\n", obligor);
    fprintf("%.4f ", intensitiesNOaccrual);
    fprintf("\nSurvival probabilities neglecting accrual for %s are:\n", obligor);
    fprintf("%.4f ", survProbs);
    fprintf("\n\n");

    stairs(0:7, [intensitiesNOaccrual; intensitiesNOaccrual(end)], 'k', 'LineWidth', 1);
    title(['CDS Bootstrap ', obligor]);
    ylabel('Intensities');
    xlabel('Years');
    hold on;

    % Compute and plot considering accrual
    flag = 2;
    [~, survProbsAccrual, intensitiesAccrual] = bootstrapCDS(dates, discounts, datesCDS, spreads, flag, R); 
    fprintf("Intensities considering accrual for %s are:\n", obligor);
    fprintf("%.4f ", intensitiesAccrual);
    fprintf("\nSurvival probabilities considering accrual for %s are:\n", obligor);
    fprintf("%.4f ", survProbsAccrual);
    fprintf("\n\n");

    stairs(0:7, [intensitiesAccrual; intensitiesAccrual(end)], 'r', 'LineWidth', 1);

    % Compute and plot JT
    flag = 3;
    [~, survProbsJT, intensitiesJT] = bootstrapCDS(dates, discounts, datesCDS, spreads, flag, R); 
    fprintf("Intensities using JT for %s are:\n", obligor);
    fprintf("%.4f ", intensitiesJT);
    fprintf("\nSurvival probabilities using JT for %s are:\n", obligor);
    fprintf("%.4f ", survProbsJT);
    fprintf("\n\n");

    stairs(0:7, [intensitiesJT; intensitiesJT(end)], 'b', 'LineWidth', 1);

    % JT check plot matches the means of intensities up to T
    Testvector = zeros(length(intensitiesAccrual), 1);
    for i = 1:length(intensitiesAccrual)
        Testvector(i) = mean(intensitiesAccrual(1:i)); 
    end
    
    stairs(0:7, [Testvector; Testvector(end)], 'g', 'LineWidth', 1);
    legend('Intensity without accrual', 'Intensity with accrual', 'Intensity J&T', 'Accrual intensity mean up to t', 'Location', 'best');
    title(['Intensities - ', obligor]);
    grid on;
    hold off;
    
    % Continuous extension of probability
    figure();
    intensitiesNOaccrual = [0; intensitiesNOaccrual];
    intensities_accrual = [0; intensitiesAccrual];
    intensities_JT = [0; intensitiesJT];
    
    datesCDS_complete = [dates(1); datesCDS]; 
    yearfracVector = [0; yearfrac(datesCDS_complete(1:end-1), datesCDS_complete(2:end), Act365)];
    
    for i = 2:length(intensitiesNOaccrual)
        dt = yearfrac(datesCDS_complete(1), datesCDS_complete(i-1), Act365); 
        f_without_accrual = @(x) exp((sum(-intensitiesNOaccrual(1:i-1).*yearfracVector(1:i-1)) - intensitiesNOaccrual(i)*(x-dt))); 
        f_accrual = @(x) exp((sum(-intensities_accrual(1:i-1).*yearfracVector(1:i-1)) - intensities_accrual(i)*(x-dt)));
        f_JT = @(x) exp((sum(-intensities_JT(1:i-1).*yearfracVector(1:i-1)) - intensities_JT(i)*(x-dt)));
        fplot(f_without_accrual, [-2+i, i-1], "k");
        hold on;
        fplot(f_accrual, [-2+i, i-1], "r");
        fplot(f_JT, [-2+i, i-1], "b");
    end
    
    grid on;
    plot(0:7, [1; survProbs], 'ok');
    plot(0:7, [1; survProbsAccrual], 'or');
    plot(0:7, [1; survProbsJT], 'ob');
    legend("Survival prob without accrual", "Survival prob with accrual", "Survival prob J&T");
    title(['Continuous Probability from Intensities - ', obligor]);
    hold off;
end


%% First to Default Pricing

% Set the dates and the CDS spreads for the two obligors useful for the
% problem
datesCDS=datesSet.swaps(1:4);
CDS_spreads_ISP=[29, 32, 35, 39]*10^(-04);
CDS_spreads_UCG=[34, 39, 45, 46]*10^(-04);
spreadsCDS=[CDS_spreads_ISP', CDS_spreads_UCG'];

% Set the expiry of the first to default (i.e the index of dates after the
% present)
expiry=5;

% Set the recovery rates
R=[0.40, 0.45];

% Set the correlation 
rho=0.2;

% Set MC confidence level
alpha=0.05;

% Set the number of MC simulations
M=100000;

% Set the seed
seed=4;

% MC simulations
[spreadmean, confidence_interval, elapsed_time]=simulateFTD(dates, discounts, datesCDS, spreadsCDS, rho, R, expiry, M, seed, alpha)

% Print the results
fprintf('Simulation Results:\n');
fprintf('Mean Spread: %.4f basis points\n', spreadmean * 10^4);
fprintf('Confidence Interval: [%.4f, %.4f] basis points\n', confidence_interval(1) * 10^4, confidence_interval(2) * 10^4);
fprintf('Elapsed Time: %.4f seconds\n', elapsed_time);

%% Vary the number of MC simulations

% Varying the number of MC simulations
M_values = 10.^[2:6];

% Initialize arrays to store results
spreadmeans = zeros(length(M_values), 1);
confidence_intervals = zeros(length(M_values), 2);
elapsed_times = zeros(length(M_values), 1);

% Loop through each value of M and perform simulations
for i = 1:length(M_values)
    M = M_values(i);
    [spreadmeans(i), confidence_intervals(i,:), elapsed_times(i)] = simulateFTD(dates, discounts, datesCDS, spreadsCDS, rho, R, expiry, M, seed, alpha);
  
end

% Create a table of results
results_table = table(M_values', ...
                      spreadmeans * 10^4, ...
                      confidence_intervals(:, 1) * 10^4, ...
                      confidence_intervals(:, 2) * 10^4, ...
                      elapsed_times, ...
                      'VariableNames', {'M', 'MeanSpread_bp', 'CI_Lower_bp', 'CI_Upper_bp', 'ElapsedTime_sec'});

% Display the table
disp(results_table);


% Plot the results
figure;
errorbar(M_values, spreadmeans * 10^4, ...
         (spreadmeans - confidence_intervals(:, 1)) * 10^4, ...
         (confidence_intervals(:, 2) - spreadmeans) * 10^4, 'o-');
set(gca, 'XScale', 'log');
xlabel('Number of Monte Carlo Simulations (M)');
ylabel('Spread FTD (bps)');
title('Spread FTD and Confidence Intervals vs. Number of MC Simulations');
grid on;

% Plot elapsed time
figure;
plot(M_values, elapsed_times, 'o-');
set(gca, 'XScale', 'log');
xlabel('Number of Monte Carlo Simulations (M)');
ylabel('Elapsed Time (seconds)');
title('Elapsed Time vs. Number of MC Simulations');
grid on;
                
%% Spread as function of rho

% Create vector of correlations
rho=linspace(-0.99, 0.99,10);

% Initialize spreads and confidence interval
spreadmean=zeros(length(rho), 1);
confidence_interval=zeros(length(rho),2);
elapsed_time=zeros(length(rho), 1);

% Number of MC simulatios
M=100000;

% Compute the mean and confidence intervals
for jj=1:length(rho)
    [spreadmean(jj), confidence_interval(jj,:), elapsed_time(jj)]=simulateFTD(dates, discounts, datesCDS, spreadsCDS, rho(jj), R, expiry, M, seed, alpha);
end


% Plot the graph as function of rho
figure;
title("CDS spread as function of rho")
plot(rho, spreadmean*10^4, 'o-', LineWidth=2);
hold on;
plot(rho, confidence_interval(:,1)*10^4, 'o-', LineWidth=2);
plot(rho, confidence_interval(:,2)*10^4, 'o-', LineWidth=2);
legend("Spread mean (bps)", "Lower CI bound for spread mean(bps)", "Upper CI bound for spread mean(bps)", Location="southwest")
xlabel("Correlation");
ylabel("CDS spread");
grid on;

% Print total elapsed time
fprintf('Total elapsed time: %.2f\n', sum(elapsed_time));
