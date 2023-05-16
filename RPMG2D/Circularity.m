function [XMIC , YMIC , RMIC , XMCC , YMCC , RMCC , CircularityFactor] = Circularity(BW,MICMethod,MCCMethod)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ;
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;

  % find centre of shape
  S = regionprops(BW,'Centroid') ; % Center must be inside the shape

  % MIC Model
  Threshold_Clustering = 5 ;
  Neighbor_Pixels = 1 ;
  Threshold_Efficiency = 0.90 ;
  Verbose = false ;
  [XMIC , YMIC , RMIC , EMIC] = MIC_Finder(X,Y,S,BW,Threshold_Clustering,Neighbor_Pixels,Threshold_Efficiency,MICMethod,Verbose) ;


  % find New Distance with New XMIC and YMIC
  for i = 1:length(X)
    DMCC(i,1) = pdist([XMIC YMIC ; X(i) Y(i)]) ;
  end

  % Smooth Distance (Optional)
  xD = (1:length(DMCC))' ;
  [yD , ~ , ~ , ~] = csaps(xD,DMCC, 0.1, xD, w=[]) ;
  DMCC = yD ;

  % Rearrange Distance and Boundary
  base_ind = find(DMCC == min(DMCC)) ;
  XX = [X(base_ind(1)+1:end) ; X(1:base_ind(1))] ;
  YY = [Y(base_ind(1)+1:end) ; Y(1:base_ind(1))] ;
  DD = [DMCC(base_ind(1)+1:end) ; DMCC(1:base_ind(1))] ;
  X = XX ; Y = YY ; DMCC = DD ;

  % MCC Model
  Threshold_Clustering = 3 ;
  Neighbor_Pixels = 1 ;
  Verbose = false ;
  [XMCC , YMCC , RMCC] = MCC_Finder(X,Y,BW,Threshold_Clustering,Neighbor_Pixels,DMCC,MCCMethod,Verbose) ;

  CircularityFactor = sqrt(RMIC/RMCC) ;

end
