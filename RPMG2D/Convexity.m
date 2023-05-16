function [X , Y , Shape_Perimeter , Convex_Perimeter , Convex_Boundary , ConvexityFactor] = Convexity(BW)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Perimeter','ConvexHull','ConvexImage') ; % Center must be inside the shape
  s = regionprops(S.ConvexImage,'Perimeter') ;
  
  Shape_Perimeter = S.Perimeter ;
  Convex_Perimeter = s.Perimeter ;
  Convex_Boundary = S.ConvexHull ;
  
  ConvexityFactor = s.Perimeter / S.Perimeter ;

end
