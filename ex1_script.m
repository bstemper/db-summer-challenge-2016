%%      Exercise (1)
%%      
%%      Computing daily VaR(99%) for stock-only portfolio.

clear;
format compact;
load('trading_days.mat');
tic;

%% Inputs

% position: holdings of assets (1x15 row vector)

% Position for Bayer Example
% position = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

% Position for exercise (i)
% position = [0,0,0,0,0,0,0,0,0,0,42,0,0,0,0];
 
% Position for exercise (ii)
 position = [41,25,43,35,45,49,58,59,41,36,43,44,45,30,37];

% Position for exercise (iii)
% position = [1,-2,3,-4,5,-6,7,-8,9,-10,11,-12,13,-14,15];


% data_file: market data CSV file to read
data_file = 'market_data.csv';

% conf_level: confidence level of VaR (0 < decimal < 1)
conf_level = 0.99;


%% Preprocessing

% reading in financial data and transposing so that rows correspond to assets
prices = transpose(csvread(data_file,1,1));

% truncating matrix 'prices' so it doesn't include VDAX and DAX anymore
prices = prices(1:15,:);
log_prices = log(prices);

% initializing daily log changes
[number_assets, number_trading_days] = size(log_prices);
daily_log_changes = zeros(number_assets, number_trading_days);

% note that changes from day 1 to day 2 are stored in index 2
for i = 2:number_trading_days
    daily_log_changes(:,i) = log_prices(:,i) - log_prices(:,i-1);
end

%% VaR computation

VaR = zeros(number_trading_days-261,1);

% Iterate through enumerated trading days, so 
%   day 1 is 2010-01-01
%   day 262 is 2011-01-03.

for T = 262:number_trading_days
    
    % No glimpse into the future, first entry (all zeros) ignored.
    truncated_data = daily_log_changes(:,2:T);
   
    % weighting more recent values for mu and sigma higher than past values
    weighting = zeros(T-1,1);
    factor = 2/(T*(T-1));
    
    for k = 1:T-1
        weighting(k,1) = k*factor;
    end
    
    %%%%%%      Weighted VAR        %%%%%%%
    
    % weighting the data with a discrete probability measure
    muW = truncated_data * weighting;
    centralized_data=truncated_data - muW * ones(1,size(truncated_data,2));
    weightingExtended=diag(weighting);
    SigmaW = centralized_data*weightingExtended*transpose(centralized_data);
   
    % Assuming joint Gaussianity of log returns, we compute the parameters
    % of the Gaussian linearized loss function.
    coefficients = (position .* transpose(prices(:,T)));
    mean_positionW = coefficients * muW;
    variance_positionW = coefficients * SigmaW * transpose(coefficients);
    
    % Using inverse Gaussian to compute quantile
    stdv_positionW = sqrt(variance_positionW);
    VaR(T-261)=norminv(1-conf_level, mean_positionW, stdv_positionW);
    
end

%% CSV file writing

% turning vector into cell, than cell concatenation
VaR = num2cell(VaR);
output = [trading_days(262:end),VaR];

% using S. Fiedlers cell2csv function to write csv file (standard csv write
% only works for numeric values)
cell2csv('More_risk_more_fun_Ex1i.csv',output);
toc;
