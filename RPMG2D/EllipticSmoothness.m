function [X , Y , EX , EY , EP , SP , EllipticSmoothnessFactor] = EllipticSmoothness(BW)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Perimeter','MajoraxisLength','MinoraxisLength','Orientation','Centroid') ; % Center must be inside the shape
  
  a = S.MajorAxisLength/2 ;
  b = S.MinorAxisLength/2 ;
  
  % Equivalent ellipse Boundary
  Theta = 0:0.01:2*pi ;

  x = a * cos(Theta) ;
  y = b * sin(Theta) ;

  EX = x * sind(90 + S.Orientation) + y * cosd(90 + S.Orientation) ;
  EY = x * cosd(90 + S.Orientation) - y * sind(90 + S.Orientation) ;

  EX = EX + S.Centroid(1) ;
  EY = EY + S.Centroid(2) ;

  % Ellipse Perimeter
  EP = pi*(a+b)*(3*(((a-b)^2)/(((a+b)^2)*(sqrt(((-3*((a-b)^2))/((a+b)^2))+4)+10)))+1) ;
  
  % Shape Perimeter
  SP = S.Perimeter ;
  
  EllipticSmoothnessFactor = SP / EP ;


end