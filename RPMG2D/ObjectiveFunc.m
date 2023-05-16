function Dcal = ObjectiveFunc(x,xdata)
    K = length(x)/3 ;
    A = x(1:K) ;
    B = x(K+1:2*K) ;
    C = x(2*K+1:3*K) ;
##    Dcal = eps*ones(size(xdata)) ;
    Dcal = zeros(size(xdata)) ;
    for i = 1:K
        Dcal = Dcal + A(i).*exp(-((xdata-B(i))/C(i)).^2) ;
    end
end
