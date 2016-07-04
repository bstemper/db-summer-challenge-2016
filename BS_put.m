function y=BS_put(S0,K,T,r,sigma)

  d1=(log(S0./K)+(r+sigma.*sigma/2).*T)./(sigma.*sqrt(T));
  d2=(log(S0./K)+(r-sigma.*sigma/2).*T)./(sigma.*sqrt(T));
  y=K.*exp(-r.*T).*normcdf(-d2)-S0.*normcdf(-d1);
 end