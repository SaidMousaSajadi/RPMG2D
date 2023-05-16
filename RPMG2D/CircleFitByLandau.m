function Para = CircleFitByLandau(XY)
x = XY(:,1) ;
y = XY(:,2) ;

n = length(x) ;
p1 = 0 ; p2 = 0 ; p3 = 0 ; p4 = 0 ; p5 = 0 ; p6 = 0 ; p7 = 0 ; p8 = 0 ; p9 = 0 ; % inital 

for i = 1:n 
    p1 = p1 + x(i) ;
    p2 = p2 + x(i)*x(i) ;
    p3 = p3 + x(i)*y(i) ;
    p4 = p4 + y(i) ;
    p5 = p5 + y(i)*y(i) ;
    p6 = p6 + x(i)*x(i)*x(i) ;
    p7 = p7 + + x(i)*y(i)*y(i) ;
    p8 = p8 + + y(i)*y(i)*y(i) ;
    p9 = p9 + + x(i)*x(i)*y(i) ;
end
a1 = 2 * (p1*p1 - n*p2) ;
b1 = 2 * (p1*p4 - n*p3) ;
a2 = b1 ;
b2 = 2 * (p4*p4 - n*p5) ;

c1 = p2*p1 - n*p6 + p1*p5 - n*p7 ;
c2 = p2*p4 - n*p8 + p4*p5 - n*p9 ;

xc = (c1*b2 - c2*b1)/(a1*b2 - a2*b1) ;
yc = (a1*c2 - a2*c1)/(a1*b2 - a2*b1) ;
r = sqrt((p2 - 2*p1*xc + n*xc*xc + p5 - 2*p4*yc + n*yc*yc)/n) ;

Para = [xc , yc , r] ;


end