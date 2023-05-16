function [X , Y , X_Convex , Y_Convex , ShapeArea , ConvexArea , SolidityFactor , FullnessFactor] = Solidity(BW)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  S = regionprops(BW,'Area','ConvexArea','ConvexHull','Solidity') ;
  
  X_Convex = S.ConvexHull(:,1) ;
  Y_Convex = S.ConvexHull(:,2) ;
  
  ShapeArea = S.Area ;
  ConvexArea = S.ConvexArea ;

  SolidityFactor = ShapeArea/ConvexArea ;
  FullnessFactor = sqrt(ShapeArea/ConvexArea) ;
end