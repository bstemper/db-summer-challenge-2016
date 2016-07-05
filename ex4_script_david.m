%%      Exercise (4)
%%      
%%      Computing daily VaR(99%) for stock-only portfolio.

clear;
clf
load('trading_days.mat');
tic();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%/ Inputs %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Inputs
% 
% position_assets = [1,-2,3,-4,5,-6,7,-8,9,-10,11,-12,13,-14,15];
% 
% position_calls = [15,-14,13,-12,11,-10,9,-8,7,-6,5,-4,3,-2,1];
% moneyness_calls = [1,0.9,0.7,1.2,2,0.5,1,0.8,1,3,1.1,1.2,0.7,1.4,1.5];
% maturity_calls = [5,6,7,8,9,0.25,0.5,0.75,1,1.25,1.5,1.75,1.75,2,2.25];
% 
% position_puts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% moneyness_puts = [1,0.9,0.7,1.2,2,0.5,1,0.8,1,3,1,1.2,0.7,1.4,1.5];
% maturity_puts = [5,6,7,8,9,0.25,0.5,0.75,1,1.25,7,1.75,1.75,2,2.25];
% 
% infos_uoc = [1,1,0.5,6];



%% Excerise 4 i

% position_assets = [1,-2,3,-4,5,-6,7,-8,9,-10,11,-12,13,-14,15];
% position_calls = [15,-14,13,-12,11,-10,9,-8,7,-6,5,-4,3,-2,1];
% moneyness_calls = [1,0.9,0.7,1.2,2,0.5,1,0.8,1,3,1.1,1.2,0.7,1.4,1.5];
% maturity_calls = [5,6,7,8,9,0.25,0.5,0.75,1,1.25,1.5,1.75,1.75,2,2.25];
% position_puts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% moneyness_puts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% maturity_puts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% infos_uoc = [1,1,0.5,6];

%Exercise 4 ii
% 
% position_assets = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% position_calls = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% moneyness_calls = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% maturity_calls = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% infos_uoc = [20,1,0.8,6];
% position_puts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% moneyness_puts = [1,0.9,0.7,1.2,2,0.5,1,0.8,1,3,1,1.2,0.7,1.4,1.5];
% maturity_puts = [5,6,7,8,9,0.25,0.5,0.75,1,1.25,7,1.75,1.75,2,2.25];

%Exercise 4 iii

position_assets = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
position_calls = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
moneyness_calls = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
maturity_calls = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
position_puts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
moneyness_puts = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
maturity_puts = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
infos_uoc = [1,1,0.7,7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Be aware! If we have no puts, uoc or calls then we have to set moneyness
% and maturity to zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data_file: market data CSV file to read
data_file = 'market_data.csv';

% conf_level: confidence level of VaR (0 < decimal < 1)
conf_level = 0.99;
learning = 262;

%Number of paths of the model we simulate
N=100;
V=N-N*conf_level;

M=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reading in financial data and transposing so that rows correspond to assets
prices = transpose(csvread(data_file,1,1));

% truncating matrix 'prices' so it doesn't include VDAX and DAX anymore
dax=prices(17,:);
vdax=0.01*prices(16,:)/260;
prices = prices(1:15,:);

% transforming prices to log prices, leaving DAX & VDAX unaffected
log_dax=log(dax);
log_prices = log(prices(1:15,1:end));

% transforming moneyness to strikes (strike * spot price on 3rd Jan 2011)
% strikes_calls = moneyness_calls .* transpose(prices(1:15,262));
% strikes_puts = moneyness_puts .* transpose(prices(1:15,262));
strikes_calls = (1./moneyness_calls) .* transpose(prices(1:15,262));
strikes_puts = moneyness_puts .* transpose(prices(1:15,262));
strike_uoc = (1/infos_uoc(2)) * dax(262);

% transforming barrier moneyness to barrier (barier * spot price on 3rd Jan 2011)
barrier_uoc = (1/infos_uoc(3))* dax(262);

% initializing daily log changes
[number_assets, number_trading_days] = size(log_prices);
daily_log_changes = zeros(number_assets, number_trading_days);
daily_log_dax_changes = zeros(number_assets, number_trading_days);
daily_call_changes = zeros(number_assets, number_trading_days);
daily_put_changes = zeros(number_assets, number_trading_days);
daily_ouc_changes = zeros(number_assets, number_trading_days);

% initializing call prices (value of the call option)
call_values = zeros(number_assets,number_trading_days);
put_values = zeros(number_assets,number_trading_days);
uoc_values = zeros(1,number_trading_days);

% note that changes from day 1 to day 2 are stored in index 2
for i = 2:number_trading_days
    daily_log_changes(:,i) = log_prices(:,i) - log_prices(:,i-1);
    daily_log_dax_changes(:,i) = log_dax(i) - log_dax(i-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% VaR computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VaR = zeros(number_trading_days-1-learning,1);
Actual_Loss = zeros(number_trading_days-1-learning,1);
counter=0;

%V=var_pointer_calibration(131,262,prices(:,1:262),daily_log_changes(:,1:262),conf_level,V,M);

% Iterate through enumerated trading days, so 

for T = learning:number_trading_days
    
        % No glimpse into the future, first entry (all zeros) ignored.
        truncated_data = daily_log_changes(:,2:T);
        truncated_data_dax = daily_log_dax_changes(2:T);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% weighting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
        % weighting more recent values for mu and sigma higher than past values
        weighting = zeros(T-1,1);
        
        % Correspond to linear weighting
        Tea = 2/(T*(T-1));
             
        % Correspond to exponential weighting
        %w = 0.968;
        %V1 = (1-w^(T-1))/(1-w);
        
        for k = 1:T-1
            % Classical weighting
            %weighting(k,1)=1/(T-1); 
            
            % linear weighting
             weighting(k,1) = k*Tea;
             
            % exponential weighting
            %weighting(k,1)=w^(T-k)/V1;
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Computation of mean (variance optimized) and covariance of changes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        mu = cov(truncated_data')*diag(sum(cov(truncated_data'),2))*truncated_data * weighting;
        centralized_data=truncated_data - mu * ones(1,size(truncated_data,2));
        weightingExtended=diag(weighting);
        Sigma = centralized_data*weightingExtended*transpose(centralized_data);
           
        % Computing annualized vol
        ann_vol = zeros(15,1);
        for i=1:15
            ann_vol(i) = sqrt(Sigma(i,i)*260);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Computation of mean of dax
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        mu_dax = cov(truncated_data_dax')*diag(sum(cov(truncated_data_dax'),2))*truncated_data_dax * weighting;
            
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Call/Put Prices (calculating call prices with explicit BS-formula)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        time_to_maturity_calls=max(maturity_calls*260+262-T,0); 
        time_to_maturity_puts=max(maturity_puts*260+262-T,0); 
        time_to_maturity_uoc=max(infos_uoc(4)*260+262-T,0);
        
        call_values(:,T)=   BS_call(prices(:,T),strikes_calls(:),time_to_maturity_calls(:),0,diag(Sigma));
        put_values(:,T)=BS_put(prices(:,T),strikes_puts(:),time_to_maturity_puts(:),0,diag(Sigma));
        uoc_values(T)=UO_call(dax(T),strike_uoc,time_to_maturity_uoc,barrier_uoc,vdax(T));
               
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Simulating the value of tommorow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A=chol(Sigma);
        Y=normrnd(0,1,number_assets,N);
        W=A'*Y;
        X=W+(diag(mu)-1/2*diag(diag(Sigma)))*ones(number_assets,N);
        Future=diag(prices(:,T))*exp(X);
        
        X_dax=vdax(T)*normrnd(0,1,1,N)+mu_dax*ones(1,N)-1/2*vdax(T)^2*ones(1,N);
        Future_dax=dax(T)*exp(X_dax);   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calculating the Loss of tommorow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Loss_assets=position_assets*(Future(:,1:N)-prices(:,T)*ones(1,N));
        Loss_calls=position_calls*(BS_call(Future,strikes_calls(:)*ones(1,N),max((time_to_maturity_calls(:)-1),0)*ones(1,N),0,(diag(Sigma)+ann_vol/260)*ones(1,N))-call_values(:,T)*ones(1,N));  
        Loss_puts=position_puts*(BS_put(Future,strikes_puts(:)*ones(1,N),max((time_to_maturity_puts(:)-1),0)*ones(1,N),0,(diag(Sigma)+ann_vol/260)*ones(1,N))-put_values(:,T)*ones(1,N));  
        Loss_uoc=infos_uoc(1)*(UO_call(Future_dax,strike_uoc,max((time_to_maturity_uoc-1),0),barrier_uoc,vdax(T))-uoc_values(T));
        Loss=Loss_assets+Loss_calls+Loss_puts+Loss_uoc;
        
        Sorted_Loss=sort(Loss);
        VaR(T-learning+1)= Sorted_Loss(V);
        if T== number_trading_days
            Actual_Loss(T-learning+1)=0;
        else
            Actual_Loss_assets=position_assets * (prices(:,T+1)-prices(:,T))
            
            Actual_Loss_calls=position_calls*(BS_call(prices(:,T+1),strikes_calls(:),max(time_to_maturity_calls(:)-1,0),0,diag(Sigma))-call_values(:,T))
            Actual_Loss_puts=position_puts*(BS_put(prices(:,T+1),strikes_puts(:),max(time_to_maturity_puts(:)-1,0),0,diag(Sigma))-put_values(:,T))
            Actual_Loss_uoc=infos_uoc(1)*(UO_call(dax(T+1),strike_uoc,max(time_to_maturity_uoc-1,0),barrier_uoc,vdax(T))-uoc_values(T))
            Actual_Loss(T-learning+1) = Actual_Loss_assets+Actual_Loss_calls+Actual_Loss_puts+Actual_Loss_uoc;
        end
        
        %Calculating the fails of our VaR on the given data set    
        if VaR(T-learning+1)>Actual_Loss(T-learning+1)
            counter=counter+1;
        end      
end

x=(1:1:number_trading_days+1-learning);
plot(x,VaR,x,Actual_Loss);
real_Var=1-counter/(number_trading_days-learning)
toc();
