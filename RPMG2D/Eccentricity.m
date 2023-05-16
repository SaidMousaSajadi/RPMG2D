function [a , b , c , Ellipse_Boundary , MajorX1 , MajorX2 , MajorY1 , MajorY2 , MinorX1 , MinorX2 , MinorY1 , MinorY2 , XCent , YCent , cx , cy , bx , by , EccentricityFactor] = Eccentricity(BW) 
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Eccentricity','MajoraxisLength','MinoraxisLength','Orientation','Centroid') ; % Center must be inside the shape
  
  # a^2 - b^2 = c^2
  a = S.MajorAxisLength/2 ;
  b = S.MinorAxisLength/2 ;
  c = sqrt(a^2 - b^2) ;
  
  Theta = 0:0.01:2*pi ;

  % step 1st :
  xelli = a * cos(Theta) ;
  yelli = b * sin(Theta) ;
  % step 2nd :
  XElli = xelli*sind(90 + S.Orientation) + yelli*cosd(90 + S.Orientation) ;
  YElli = xelli*cosd(90 + S.Orientation) - yelli*sind(90 + S.Orientation) ;
  % step 3th : Perimeter of equivalent ellipse
  XElli = XElli + S.Centroid(1) ;
  YElli = YElli + S.Centroid(2) ;
  
  % Equivalent ellipse Boundary
  Ellipse_Boundary = [XElli ; YElli] ;
  
  XCent = S.Centroid(1) ;
  YCent = S.Centroid(2) ;
  
  MajorX1 = S.Centroid(1) + a*cosd(S.Orientation) ;
  MajorY1 = S.Centroid(2) - a*sind(S.Orientation) ;
  MajorX2 = S.Centroid(1) - a*cosd(S.Orientation) ;
  MajorY2 = S.Centroid(2) + a*sind(S.Orientation) ;
  
  MinorX1 = S.Centroid(1) + b*cosd(90+S.Orientation) ;
  MinorY1 = S.Centroid(2) - b*sind(90+S.Orientation) ;
  MinorX2 = S.Centroid(1) - b*cosd(90+S.Orientation) ;
  MinorY2 = S.Centroid(2) + b*sind(90+S.Orientation) ;

  % On Major Diameter
  cx = XCent + c * cosd(S.Orientation) ;
  cy = YCent - c * sind(S.Orientation) ;
  
  % On Minor Diameter
  bx = XCent + b * cosd(90+S.Orientation) ;
  by = YCent - b * sind(90+S.Orientation) ;
  
  EccentricityFactor = c/a ; % S_Eccentricity = S.Eccentricity ;
##  



end