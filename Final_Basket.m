clear; clc;

% 1. Retrieve the data for your empirical analysis (e.g., from Datastream, 
% Bloomberg, and Reuters). Use monthly data as base case.

% Import bloomberg data and split up in needed formats
[ Data, Headertext ] = xlsread('C:\Users\Cedric\Desktop\Uni\19_FS\Financial Modeling - Asset Allocation\Homework\AssetAllocation.xlsx')
% Define date vector for plots later on
date = datenum( Headertext( 3:end , 1 ), 'dd.mm.yyyy' );
% Headertexts as reference
Headertext = Headertext( 2, 2:end );
% Initial Investment of 1 Mio for evaluating absolute values
Init_Inv = 1000000;
% Define Returns and stock prices
Libor_1M = Data( 1:(end - 1) , end ) / 1200; % Libor rate
Return_Libor_1M = log ( 1 + Libor_1M ); % Libor return
T_Bill_1M = Data( 1:(end - 1) , end - 1 )/1200; % T-Bill rate
Return_T_Bill_1M = log (1 + T_Bill_1M ); % T-Bill return
Market = Data( : , end - 2 ); % S&P 500 index values
Data = Data ( : , 1 : (end - 3) ); % Stock prices


% 2. Find a portfolio construction mechanism and implement it (e.g., 
% Mean-Variance Portfolio optimization, equally weighted portfolio, etc.). 
% Comment on the weight restrictions used in your analysis. 
% Deliver the precise formulas used to come up with the portfolio weights 
% and elaborate on rebalancing frequency.

% Create n_shares, which is the matrix that stores the number of each share
% per month (i.e. monthly rebalancing) for an equally-weighted portfolio
n_assets = size( Data, 2 );
n_shares = ( n_assets ) ./ Data;



% 3. Determine the performance of the basket given the chosen investment 
% universe and portfolio construction methodology in the last 5-10 years.

% 7. Additionally, the analysis should consider exchange rate risk and 
% transaction costs where applicable. Ensure that exchange rate effects 
% are properly accounted for in your performance calculations and describe 
% concisely how this is achieved.

% Returns w/o including Transaction Costs, i.e. "Pure returns"
Pure_Returns = diff( log( Data ) );
Pure_Monthly_Performance = (sum( Pure_Returns' ) / n_assets)';
Pure_Cumm_Performance = cumprod( 1 + Pure_Monthly_Performance );
Pure_Abs_Performance = Pure_Cumm_Performance * Init_Inv ;


%Calculate Trading costs (i.e. Spread, Brokerage fee and Price impact)
Perc_Trading_costs= 0.0042;
Traded_Shares = abs( diff( n_shares ) );
Value_Traded_Shares = ( Traded_Shares .* Data( 2:end , : ) );
Portfolio_Value = Data .* n_shares;

% Turnover is defined the OGAW directive as sum((Value Shares Sold, 
% Value Shares Bought)/ Portfolio Value.
Rel_Portf_Turnover = sum( Value_Traded_Shares' )' ./ sum( Portfolio_Value( 2:end , : )' )';
Trading_Costs = Rel_Portf_Turnover  * Perc_Trading_costs;


% Return including Trading Costs
Adj_Monthly_Performance = (( sum( Pure_Returns' - Trading_Costs' ) ) / n_assets )';
Adj_Cumm_Performance = cumprod( 1 + Adj_Monthly_Performance );
Adj_Abs_Performance = Adj_Cumm_Performance * Init_Inv;

%Monthly transaction costs (or fees) of 0.15%
Transaction_Costs = 0.0015;
% Return including Transaction Costs
Fully_Adj_Monthly_Performance = ( sum(( Pure_Returns' - Transaction_Costs )) / n_assets )';
Fully_Adj_Cumm_Performance = cumprod( 1 + Fully_Adj_Monthly_Performance );
Fully_Adj_Abs_Performance = Fully_Adj_Cumm_Performance * Init_Inv;

% Plot the Cummulative returns
Performance = [ [ 1 ; Pure_Cumm_Performance ] [ 1 ; Adj_Cumm_Performance ]  ...
        [ 1 ; Fully_Adj_Cumm_Performance ]]
plot( date, Performance );
legend( 'Performance', 'Performance less Trading Costs', ...
    'Performance less Fees');
xlabel( 'Year' );
ylabel( 'Performance' );
xlim([ date(1) date(end) ])
datetick( 'x' , 11 , 'keeplimits' );

% Plot hist basket performance against S&P
Returns_Market = diff( log( Market ) );
Cumm_Market_Perf = cumprod( 1 + Returns_Market );
Performance = [ [ 1 ; Pure_Cumm_Performance ] [ 1 ; Cumm_Market_Perf ] ]
plot( date, Performance );
legend( 'Equity Basket', 'S&P 500');
xlabel( 'Year' );
ylabel( 'Performance' );
xlim([ date(1) date(end) ])
datetick( 'x' , 11 , 'keeplimits' );



% Matrix for absolute return with Init_Inv = 1 Mio.
Abs_Performance =[ [ Init_Inv ; Pure_Abs_Performance] [Init_Inv ; ...
    Adj_Abs_Performance] [Init_Inv ; Fully_Adj_Abs_Performance ]];


% 4. Calculate the basket risk characteristics, e.g. Value-at-Risk, 
% shortfall measures, volatility...

% Continue with Pure Performance

% Historic Variance-Covariance Matrix
Covariance_Matrix_Historic = cov( Pure_Returns );
Correlation_Matrix_Historic = corr( Pure_Returns );

% Portfolio weights
Portfolio_Weights = ones( 1 , size( Pure_Returns, 2 ) ) / size( Pure_Returns, 2 );

% Portfolio Variance / Standard Deviation
Portfolio_Variance = ( Portfolio_Weights * Covariance_Matrix_Historic * Portfolio_Weights' ) *12;
Portfolio_Std_Dev = sqrt( Portfolio_Variance );

%  Annualized VAR Variance-Covariance Methodology for 1% and 5% perc
VAR_Covar = [ ( mean( Pure_Monthly_Performance ) - Portfolio_Std_Dev / sqrt(12)...
    * 1.645  ) ( mean( Pure_Monthly_Performance ) - Portfolio_Std_Dev / sqrt(12) * 2.326 ) ];

% Monthly VAR Historical Methodology for 1% and 5% perc
Sorted_Month_Perf = sort( Pure_Monthly_Performance );
VAR_Hist = prctile( Pure_Monthly_Performance , [5 1] );
% or alternatively without linear interpolation
Position_95 = round( size( Pure_Returns, 1 ) / 20 );
Position_99 = round( size( Pure_Returns, 1 ) / 100 );
Sorted_Month_Perf = sort( Pure_Monthly_Performance );
VAR_Hist = [ (Sorted_Month_Perf( Position_95, : ) ) ...
        (Sorted_Month_Perf( Position_99, : ) ) ];
    
% Monthly VAR Simulation Methdology for 1% and 5% perc

Random_Numbers = normrnd( 0 , 1 , 100000 , 1 );
Portfolio_Returns_Simulated = (mean( Pure_Monthly_Performance ) + ...
    Portfolio_Std_Dev * Random_Numbers / sqrt(12)) 
VAR_Simulated = prctile( Portfolio_Returns_Simulated , [ 5 1 ] )
    
% Expected shortfall
Sorted_Month_Perf = sort( Pure_Monthly_Performance );
Expected_Shortfall = [ (mean( Sorted_Month_Perf( 1 : Position_95, : )) ) ...
    (mean( Sorted_Month_Perf( 1 : Position_99, : )) ) ];

% Percentage contribution to risk
Portfolio_Returns = Pure_Returns * Portfolio_Weights';
Covariance_Assets_Portfolio_help = cov( [ Portfolio_Returns Pure_Returns ] );
Covariance_Assets_Portfolio = Covariance_Assets_Portfolio_help( 1, 2 : ...
    size( Covariance_Assets_Portfolio_help , 2 ) );   
Risk_Contribution = Portfolio_Weights .* Covariance_Assets_Portfolio / Portfolio_Variance;

% Create pie chart with percentage contribution to risk
Rel_Risk_Contribution = Risk_Contribution / sum( Risk_Contribution );
% needed to create labels with 1 decimal
percent = Rel_Risk_Contribution' * 100;
labels = cellstr( num2str( percent, '%.1f%%' ) );
% Create actual pie chart
pie( Rel_Risk_Contribution , labels )
legend('JNJ', 'PFE', 'MRK', 'ABT', 'GILD', 'LLY', 'AMGN', 'BMY', 'AGN', ...
'CELG', 'BIIB', 'MYL', 'BAX', 'REGN', 'HLF', 'ALXN', 'ILMN', 'UTHR');

% Maximum Drawdown
Pure_Cumm_Performance_Ext = [ 1 ; Pure_Cumm_Performance ]
% Initialize Matrix
Maximum_Drawdown = zeros( size( Pure_Cumm_Performance_Ext , 1 ), 1 )
Drawdown_help = zeros(size( Pure_Cumm_Performance_Ext , 1) , 1 )
% Calculate Drawdown for first iteration for the loop to work aftwards
Drawdown = max( Pure_Cumm_Performance_Ext( 1, 1 ) - ...
    min( Pure_Cumm_Performance_Ext( 1 : end, 1 ) ), 0 ) ...
    / Pure_Cumm_Performance_Ext( 1, 1 )
Drawdown_help(1) =  Drawdown; 
Maximum_Drawdown(1,1) = max( Maximum_Drawdown(1,1) , Drawdown );

for i = 2 : size( Data, 1 )
    Drawdown = max( Pure_Cumm_Performance_Ext( i, 1 ) - ...
           min( Pure_Cumm_Performance_Ext( i : end, 1 ) ), 0 ) ...
           / Pure_Cumm_Performance_Ext( i, 1 )
    Drawdown_help(i) =  Drawdown  
    Maximum_Drawdown(i,1) = max( Maximum_Drawdown(i-1,1), Drawdown )
end

% Plot Max DD as well as DD
DD = [ Maximum_Drawdown Drawdown_help]
plot( date, DD );
legend( 'Maximum Drawdown' , 'Drawdown' );
xlabel( 'Year' );
ylabel( 'Drawdown in %' );
xlim([date(1) date(end)])
datetick( 'x' , 11 , 'keeplimits' );


% 5. Calculate performance measures for the basket. Among others: Mean 
% returns, return skewness, return kurtosis, Sharpe Ratio, Treynor Ratio...

% Portfolio mean / median / skewness / kurtosis return
Portfolio_Mean = mean( Pure_Monthly_Performance );
Portfolio_Median = median( Pure_Monthly_Performance );
Portfolio_Skewness = skewness( Pure_Monthly_Performance );
Portfolio_Kurtosis = kurtosis( Pure_Monthly_Performance );

% Stock mean / median / skewness / kurtosis returns
Stock_Mean = mean ( Pure_Returns );
Stock_Median = median ( Pure_Returns );
Stock_Skewness = skewness( Pure_Returns );
Stock_Kurtosis = kurtosis( Pure_Returns );

% Market mean / median / skewness / kurtosis returns
Returns_Market = diff( log( Market ) );                  
Market_Mean = mean ( Returns_Market );
Market_Median = median ( Returns_Market );
Market_Skewness = skewness( Returns_Market );
Market_Kurtosis = kurtosis( Returns_Market );


% Sharpe Ratio Portfolio
Sharpe_Ratio_Portfolio = ( mean( ( Pure_Monthly_Performance ) - Return_Libor_1M ) * 12) ...
        / Portfolio_Std_Dev;
Sharpe_Ratio_Market = ( mean( ( Returns_Market ) - Return_Libor_1M ) * 12 )...
        / Portfolio_Std_Dev;


% Jensen's Alpha
% Calculating monthly market returns and then market excess returns
Excess_Returns_Market = Returns_Market - Return_Libor_1M; 

% Now estimate alpha and beta of portfolio
Coeffs = regress( ( Pure_Monthly_Performance ), ...
[ ones( size( Excess_Returns_Market , 1 ) , 1 ) Excess_Returns_Market ])
Jensens_Alpha = Coeffs( 1 , : ) * 12 ;
Beta = Coeffs( 2 , : );

% Treynor Ratio
Treynor_Ratio = (mean( ( Pure_Monthly_Performance ) - Return_Libor_1M ) * 12) / Beta;

% Tracking Error
% First get R-Squared
mdl = fitlm( Excess_Returns_Market , ( Pure_Monthly_Performance ) );
R_squared= mdl.Rsquared.Ordinary;
% Then compute Tracking Error
Tracking_Error = Portfolio_Std_Dev * sqrt( 1 - R_squared );

% Appraisal Ratio or Treynor/Black Ratio
Appraisal_Ratio = Jensens_Alpha / Tracking_Error;


% 6. Determine the factor exposure of the basket using the Fama?French risk
% factors. These factors are available for download at 
% "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html”
% in the “Historical Benchmark Returns (Downloadable Files)” section. 
% Choose “Fama/French Benchmark Factors.”

% Read-in modified (deleted all months that are not required and changed 
% format to xlsx) FF_Benchmark file
[ FF_Factors, Headertext_2 ] = xlsread('C:\Users\Cedric\Desktop\Uni\19_FS\Financial Modeling - Asset Allocation\Homework\FF_Bench_Factors_Month_Mod.xlsx')
% Delete first column as it is a date and no F-F factor
FF_Factors = FF_Factors( 1 : end , 2 : end ) / 100


% Calculate factor exposures over time with 3 year rolling-window
Fund_Exposure_Results = zeros( size( Pure_Monthly_Performance , 1 ) -35 , 4);
for i = 36 : size( Pure_Monthly_Performance, 1 )
    Fund_Returns_help = Pure_Monthly_Performance(  i-35 : i, 1 ) ...
        - Return_Libor_1M( i-35 : i , 1 );
    Factors_help = FF_Factors( i-35 : i , : );
    b = regress( Fund_Returns_help , ... 
        [ ones( size( Fund_Returns_help , 1 ), 1 ) Factors_help ] );
    Fund_Exposure_Results( i-35 , : ) = b';
end

% Illustrate style drift using a graph
date2 = date( 37 : end , : )
plot( date2 , Fund_Exposure_Results );
legend( 'Intercept' , 'Market' , 'SMB' , 'HML' );
xlabel( 'Year' );
ylabel( 'Parameter Estimate' );
xlim( [ date2(1) date2(end) ] )
datetick( 'x' , 11 , 'keeplimits' );

% Average Factor Exposure over observed time frame
Average_Fund_Exposure = mean( Fund_Exposure_Results);





