function [X , Y , Perimeter , Area , CompactnessFactor] = Compactness(BW) 
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Centroid',"Perimeter","Area") ; % Center must be inside the shape
  
  Area = S.Area ;
  Perimeter = S.Perimeter ;
  CompactnessFactor = (Perimeter.^2)/Area ;  
end