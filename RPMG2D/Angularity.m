function [XMIC , YMIC , RMIC , EMIC , XFar , YFar , MaxDist , XCircMIC , YCircMIC , Alpha , AngularityFactor] = Angularity(BW,Method)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ;
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;

  % find centre of shape
  S = regionprops(BW,'Centroid') ; % Center must be inside the shape

  % Check Centre is inside the shape?
  % figure ; imshow(Binary_Shape) ; hold on ; plot(S.Centroid(1) , S.Centroid(2) , "*b")
  Threshold_Clustering = 5 ;
  Neighbor_Pixels = 1 ;
  Threshold_Efficiency = 0.90 ;
  Verbose = false ;
  [XMIC , YMIC , RMIC , EMIC] = MIC_Finder(X,Y,S,BW,Threshold_Clustering,Neighbor_Pixels,Threshold_Efficiency,Method,Verbose) ;



  % Distance between MIC center and boundary
  for i = 1:length(X)
    Dis(i) = pdist([XMIC YMIC ; X(i) Y(i)]) ;
  end
  MaxDist = max(Dis) ;

  % Find farthest point
  Alpha = 2*asind(RMIC/MaxDist) ; % degree
  MaxInd = find(Dis == MaxDist) ;
  XFar = X(MaxInd(1)) ;
  YFar = Y(MaxInd(1)) ;

  % Distance tangent point of circle to farthest point from MIC
  Side = sqrt(MaxDist^2 - RMIC^2) ;
  theta = 0:0.01:2*pi ;
  xCircMIC = XMIC + RMIC*cos(theta) ;
  yCircMIC = YMIC + RMIC*sin(theta) ;

  % Find tangent point of circle to farthest point from MIC
  for i = 1:length(xCircMIC)
    Sider(i) = pdist([XFar YFar ; xCircMIC(i) yCircMIC(i)]) ;
  end
  DSide = abs(Sider - Side) ;
  PointInd = find(DSide == min(DSide)) ;
  XCircMIC = xCircMIC(PointInd(1)) ;
  YCircMIC = yCircMIC(PointInd(1)) ;

  % Calculate Angularity
  AngularityFactor = (180 - Alpha)*(MaxDist/RMIC) ;

end
