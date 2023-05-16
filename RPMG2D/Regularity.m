function [X , Y , X_Convex , Y_Convex , Shape_Perimeter , Convex_Perimeter , RegularityFactor] = Regularity(BW)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Perimeter','ConvexHull','ConvexImage') ; % Shape % Center must be inside the shape
  s = regionprops(S.ConvexImage,'Perimeter') ; % Convex
  
  X_Convex = S.ConvexHull(:,1) ;
  Y_Convex = S.ConvexHull(:,2) ;
  
  Shape_Perimeter = S.Perimeter ;
  Convex_Perimeter = s.Perimeter ;
  
  RegularityFactor = log(Shape_Perimeter/(Shape_Perimeter-Convex_Perimeter)) ;

end
