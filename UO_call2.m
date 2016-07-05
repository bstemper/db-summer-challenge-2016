function y=UO_call2(S,K,T,sigma,B)

    if S>=B
       y=0;
      else
      mu=-(sigma*sigma)/2;
      a=mu/(sigma*sigma);
      b=mu/(sigma*sigma);
      lambda=1+(mu/(sigma*sigma));
      z=(log(B/S))/(sigma*sqrt(T))+b*sigma*sqrt(T);
      
      
      x1=log(S/K)/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
      x2=log(S/B)/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
      y1=log((B*B)/(S*K))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
      y2=log(B/(S))/(sigma*sqrt(T))+lambda*sigma*sqrt(T);
      z=log(B/S)/(sigma*sqrt(T))+b*sigma*sqrt(T);
      eta=-1;
      phi=1;
  
      I1=phi*S*normcdf(phi*x1)-phi*K*normcdf(phi*x1-phi*sigma*sqrt(T));
      I2=phi*S*normcdf(phi*x2)-phi*K*normcdf(phi*x2-phi*sigma*sqrt(T));
      I3=phi*S*((B/S)^(2*lambda))*normcdf(eta*y1)-phi*K*((S/B)^(2*lambda-2))*normcdf(eta*y1-eta*sigma*sqrt(T));
      I4=phi*S*((B/S)^(2*lambda))*normcdf(eta*y2)-phi*K*((S/B)^(2*lambda-2))*normcdf(eta*y2-eta*sigma*sqrt(T));
      I6=((B/S)^(a+b))*normcdf(eta*z)-((B/S)^(a-b))*normcdf(eta*z-2*eta*b*sigma*sqrt(T));
  
   end
    if K>B
        y=max(I6,0);
    else
        y=max(I1-I2+I3-I4+I6,0);
    end
 end