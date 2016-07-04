function y=UO_call(S0,K,T,sigma,B)
    x1=log(S0./K)./(sigma.*sqrt(T))+0.5*sigma.*sqrt(T);
    x2=log(S0/B)/(sigma*sqrt(T))+0.5*sigma*sqrt(T);
    y1=log(B.*B./(S0.*K))./(sigma.*sqrt(T))+0.5*sigma.*sqrt(T);
    y2=log(B./(S0))./(sigma.*sqrt(T))+0.5*sigma.*sqrt(T);
    z=log(B./S0)./(sigma.*sqrt(T))+0.5*sigma.*sqrt(T);
    eta=-1;
    phi=1;
    
    A=phi*S0.*normcdf(phi*x1)-phi*K.*normcdf(phi*x1-phi*sigma.*sqrt(T));
    B=phi*S0.*normcdf(phi*x2)-phi*K.*normcdf(phi*x2-phi*sigma.*sqrt(T));
    C=phi*S0.*(B./S0).*normcdf(eta*y1)-phi*K.*(S0./B).*normcdf(eta*y1-eta*sigma.*sqrt(T));
    D=phi*S0.*(B./S0).*normcdf(eta*y2)-phi*K.*(S0./B).*normcdf(eta*y2-eta*sigma.*sqrt(T));
    F=normcdf(eta*z)+(S0./B).*normcdf(eta*z-eta*sigma.*sqrt(T));
   
    if K>B
        y=F;
    else
        y=A-B+C-D+F;
    end
 end