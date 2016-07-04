function y=BS_call(S0,K,T,r,sigma)

  d1=(log(S0./K)+(r+sigma.*sigma/2).*T)./(sigma.*sqrt(T));
  d2=(log(S0./K)+(r-sigma.*sigma/2).*T)./(sigma.*sqrt(T));
  y=S0.*normcdf(d1)-K.*exp(-r.*T).*normcdf(d2);
 end