%%      Exercise (2)
%%      
%%      Computing daily VaR(99%) for portfolio including stocks and 
%%      European call options.

clear;
format compact;

% Loads all trading days from 2010 onwards into a cell to be used later.
load('trading_days.mat');

tic;

%% Inputs
position_shares = [1,-2,3,-4,5,-6,7,-8,9,-10,11,-12,13,-14,15];
position_calls = [15,-14,13,-12,11,-10,9,-8,7,-6,5,-4,3,-2,1];
moneyness = [1,0.9,0.7,1.2,2,0.5,1,0.8,1,3,1.1,1.2,0.7,1.4,1.5];
maturity = [5,6,7,8,9,0.25,0.5,0.75,1,1.25,1.5,1.75,1.75,2,2.25];


% data_file: market data CSV file to read
data_file = 'market_data.csv';

% conf_level: confidence level of VaR (0 < decimal < 1)
conf_level = 0.99;

%% Preprocessing & Initialization

% reading in financial data and transposing so that rows correspond to assets
prices = transpose(csvread(data_file,1,1));

% transforming prices to log prices, leaving DAX & VDAX unaffected
log_prices = log(prices(1:15,1:end));

% transforming moneyness to strikes (strike * spot price on 3rd Jan 2011)
strikes = moneyness .* transpose(prices(1:15,262));

% initializing daily log changes
[number_assets, number_trading_days] = size(log_prices);
daily_log_changes = zeros(number_assets, number_trading_days);

% initializing Greeks
delta = zeros(1,number_assets);
theta = zeros(1,number_assets);
vega = zeros(1,number_assets);

% note that changes from day 1 to day 2 are stored in index 2
for i = 2:number_trading_days
    daily_log_changes(:,i) = log_prices(:,i) - log_prices(:,i-1);
end

% Initializing output vector
VaR = zeros(number_trading_days-261,1);


%% VaR computation

for T = 262:number_trading_days

%     %% Computing weighted variance of plain prices
%     
%     % No glimpse into the future. truncated_prices = prices(1:15,1:T);
%     
%     % Computing weights used in next section. weighting = zeros(T,1);
%     factor = 2/(T*(T+1));
%     
%     for k=1:T
%         weighting(k,1) = k*factor;
%     end
%     
%     % Computing mean and variance of reweighted plain prices
%     w_prices_mean = truncated_prices * weighting; centralized_data =
%     truncated_prices - w_prices_mean * ones(1,size(truncated_prices,2));
%     weight_diag = diag(weighting); w_prices_Sigma =
%     centralized_data*weight_diag*transpose(centralized_data);
    
    %% Computing weighted mean and variance of Gaussian joint log returns.
    
    % No glimpse into the future, first entry (all zeros) ignored.
    truncated_log_returns = daily_log_changes(:,2:T);
   
    % Computing the weights used in next block.
    weighting = zeros(T-1,1);
    factor = 2/(T*(T-1));
    
    for k = 1:T-1
        weighting(k,1) = k*factor; 
    end
    
    % Computing mean and variance of reweighted log returns
    w_log_return_mean = truncated_log_returns * weighting;
    centralized_data = truncated_log_returns - w_log_return_mean * ones(1,size(truncated_log_returns,2));
    weight_diag = diag(weighting);
    w_log_return_Sigma = centralized_data*weight_diag*transpose(centralized_data);
    
    % Computing annualized vol
    ann_vol = zeros(15,1);
    for i=1:15
        ann_vol(i) = sqrt(w_log_return_Sigma(i,i)*260);
    end
        
    
    %% Preliminaries for computation of position parameters
    
    for i=1:15
        
        % Options past their maturities are kicked out
        time_to_maturity = (maturity(i) * 260 - (T-261))/260; % in years
        
        if  time_to_maturity <= 0
            
            position_calls(i) = 0;
            delta(i) = 0;
            vega(i) = 0;
            theta(i) = 0;
           
        else
            
            % Computing greeks and storing coefficients
            d1 = (log(prices(i,T)/strikes(i)) + 1/2 * ann_vol(i)^2 * ...
            time_to_maturity)/(ann_vol(i) * sqrt(time_to_maturity));
            d2 = d1 - ann_vol(i) * sqrt(time_to_maturity);
        
            delta(i) = normcdf(d1);
            vega(i) = prices(i,T) * normpdf(d1) * sqrt(time_to_maturity);
            theta(i) =  -(prices(i,T) * normpdf(d1) * ann_vol(i))/ ...
                        (2 * sqrt(time_to_maturity));    
        
        end
    end
        
    %% Computation of parameters of overall position
    
    % Calculating coefficients
    delta_coeff = (position_shares + position_calls .* delta) .* transpose(prices(1:15,T));
    theta_coeff = position_calls .* theta;
    vega_coeff = position_calls .* vega;  
    
    % Parameters of total position
    total_position_mean =   delta_coeff * w_log_return_mean + ...
                            theta_coeff * repmat(1/260,15,1) + ...
                            vega_coeff * ann_vol/260;
                        
    total_position_var = delta_coeff * w_log_return_Sigma * transpose(delta_coeff);
    
    total_position_stdv = sqrt(total_position_var);
  
    % Computing and writing VaR
    VaR(T-261,1) = norminv(1-conf_level, total_position_mean, total_position_stdv);    
     
end

%% CSV file writing

% turning vector into cell, than cell concatenation
VaR = num2cell(VaR);
output = [trading_days(262:end),VaR];

% using S. Fiedlers cell2csv function to write csv file (standard csv write
% only works for numeric values)
cell2csv('More_risk_more_fun_Ex2i.csv',output);
toc;

