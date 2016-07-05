function y=var_pointer_calibration(learning,testing,prices,daily_log_changes,alpha,V,N)
    
    position_assets=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
    VaR = zeros(testing-learning,1);
    Actual_Loss = zeros(testing-learning,1);
    counter=0;
    
    for T = learning:testing
    
        % No glimpse into the future, first entry (all zeros) ignored.
        truncated_data = daily_log_changes(:,2:T);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% weighting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
        % weighting more recent values for mu and sigma higher than past values
        weighting = zeros(T-1,1);
        
       for k = 1:T-1
            % Classical weighting
            weighting(k,1)=1/(T-1); 
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Computation of mean (variance optimized) and covariance of changes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        mu = cov(truncated_data')*diag(sum(cov(truncated_data'),2))*truncated_data * weighting;
        centralized_data=truncated_data - mu * ones(1,size(truncated_data,2));
        weightingExtended=diag(weighting);
        Sigma = centralized_data*weightingExtended*transpose(centralized_data);
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Simulating the value of tommorow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        A=chol(Sigma);
        Y=normrnd(0,1,15,N);
        W=A'*Y;
        X=W+(diag(mu)-1/2*diag(diag(Sigma)))*ones(15,N);
        Future=diag(prices(:,T))*exp(X);
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Calculating the Loss of tommorow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Loss_assets=position_assets*(Future(:,1:N)-prices(:,T)*ones(1,N));
        Sorted_Loss=sort(Loss_assets);
        VaR(T-learning+1)= Sorted_Loss(V);
        
        if T== testing
            Actual_Loss(T-learning+1)=0;
        else
            Actual_Loss(T-learning+1)=position_assets * (prices(:,T+1)-prices(:,T));            
        end
        
        %Calculating the fails of our VaR on the given data set    
        if VaR(T-learning+1)>Actual_Loss(T-learning+1)
            counter=counter+1;
        end      
    end
    real_Var=1-counter/(testing-learning);
    if real_Var<alpha
            y=var_pointer_calibration(learning,testing,prices,daily_log_changes,alpha,V-1,N);
    else
        y=V
    end
    
    
    
 end