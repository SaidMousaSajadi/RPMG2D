function [BaseLoc , MaxLoc , MaxDiff , MinLoc , MinDiff] = findextremas(D,varargin)
  if nargin == 1
    Var = 1 ;
  elseif nargin == 2
    if lower(varargin{1}) == 'correction'
      Var = diff(linspace(1,length(D),length(D))') ;
    else
      Var = 1 ;
    end
  else
    Var = 1 ;
  end

  % To Find Max / Min Locations
  Sign = sign(diff(D)./Var) ;
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
  MaxLoc = find((sign(DSign)<0)==1) ; MaxLoc(end) = length(D) ;
  MinLoc = find((sign(DSign)>0)==1) ;

  % To Find Extresmas
  Diff = -diff(D)./Var ;
  Diff2 = -diff(Diff)./Var ;

  [~,BaseLoc] = findpeaks(Diff,"DoubleSided") ;
  BaseLoc = sort(BaseLoc) ;

  MaxDiff2Loc = [] ;
  for i = 1:size(Diff2,1)-1
    if (Diff2(i) < 0) && (Diff2(i+1) > 0)
      MaxDiff2Loc = [MaxDiff2Loc ; i+1] ;
    end
  end
  MinDiff2Loc = [] ;
  for i = 1:size(Diff2,1)-1
    if (Diff2(i) > 0) && (Diff2(i+1) < 0)
      MinDiff2Loc = [MinDiff2Loc ; i+1] ;
    end
  end

  % Maximum Extremas
  MaxDiff = sort(unique(setdiff(MaxDiff2Loc,BaseLoc))) ;
  % Minmum Extremas
  MinDiff = sort(unique(setdiff(MinDiff2Loc,BaseLoc))) ;


  end
