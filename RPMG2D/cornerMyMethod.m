function Out = cornerMyMethod(bw)
  % find boundary
  [Bound,~,~] = bwboundaries(bw) ;
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;

  S = regionprops(bw) ;

  for i = 1:length(X)
    D(i,1) = pdist([S(1).Centroid(1) S(1).Centroid(2) ; X(i) Y(i)]) ;
  end

  %% Smooth distance (Optional)
  xD = (1:length(D))' ;
  [yD , ~ , ~ , ~] = csaps(xD,D, 0.95, xD, w=[]) ;
  D = yD ;

  base_ind = find(D==min(D)) ;
  XX = [X(base_ind(1)+1:end) ; X(1:base_ind(1))] ;
  YY = [Y(base_ind(1)+1:end) ; Y(1:base_ind(1))] ;
  DD = [D(base_ind(1)+1:end) ; D(1:base_ind(1))] ;
  X = XX ; Y = YY ; D = DD ;

  [MaxLoc , MinLoc] = FindExtramas_Loc(D) ;
  [BaseLoc , ~ , MaxDiff , ~ , MinDiff] = findextremas(D) ;

##  figure ; plot(D) ; hold on ;
##  plot(MaxLoc,D(MaxLoc),'*r') ;
##  plot(MaxDiff,D(MaxDiff),'*g') ;
##  plot(MinDiff,D(MinDiff),'*b') ;
##  plot(BaseLoc,D(BaseLoc),'*k') ;

  LOC = sort(unique([MaxLoc ; MaxDiff ; MinDiff])) ;
  XCON = X(LOC) ;
  YCON = Y(LOC) ;

  Out = [XCON , YCON] ;
##  hold on ; plot(XCON,YCON,'*y')

  function [MaxLoc , MinLoc] = FindExtramas_Loc(D)
    % To Find Max / Min Locations
    Sign = sign(diff(D)./1) ;
    if isrow(Sign)
      Sign = [1 Sign] ;
    else
      Sign = [1 ; Sign] ;
    end
    Sign(end) = -1 ;

    DSign = diff(Sign) ;
    if isrow(DSign)
      DSign = [DSign 0] ;
    else
      DSign = [DSign ; 0] ;
    end

    % Maximum and Minimum
    MaxLoc = find((sign(DSign)<0)==1) ;
    MinLoc = find((sign(DSign)>0)==1) ;
  end
end
