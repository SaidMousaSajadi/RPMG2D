function [Major_X1 , Major_X2 , Major_Y1 , Major_Y2 , Minor_X1 , Minor_Y1 , Minor_X2 , Minor_Y2 , Major , Minor, AspectRatioFactor] = AspectRatio(BW,varargin)

  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ;
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;

  % find centre of shape
  S = regionprops(BW,'Centroid',"Orientation") ; % Center must be inside the shape

  % Distance matrix
  Coord = [X' ; Y'] ;
  DistMatrix = squareform(pdist(Coord')) ;
  Sorted = unique(sort(DistMatrix(:))) ;

  % Check "Major" In Shape Or Not
  switch nargin
    case 2
      if ischar(varargin{1})
        switch lower(varargin{1})
          case 'inside' % Check major axis that is inside shape
            [r,c] = MajorAxisinShape(X,Y,BW,DistMatrix,Sorted,0.98) ;
          otherwise
            [r,c] = find(DistMatrix == Sorted(end)) ;
          end
      else
        [r,c] = find(DistMatrix == Sorted(end)) ;
      end
    case 3
      if ischar(varargin{1}) & isnumeric(varargin{2})
        switch lower(varargin{1})
          case 'inside' % Check major axis that is inside shape
            [r,c] = MajorAxisinShape(X,Y,BW,DistMatrix,Sorted,varargin{2}) ;
          otherwise
            [r,c] = find(DistMatrix == Sorted(end)) ;
        end
      else
        [r,c] = find(DistMatrix == Sorted(end)) ;
      end
    otherwise
      [r,c] = find(DistMatrix == Sorted(end)) ;
  end


  % Coordinates Major
  Major_X1 = X(r(1)) ;
  Major_X2 = X(c(1)) ;
  Major_Y1 = Y(r(1)) ;
  Major_Y2 = Y(c(1)) ;

  Major_Slope = (Major_Y2 - Major_Y1) / (Major_X2 - Major_X1) ;
  Minor_Slope = -1/Major_Slope ;

  % Coordinates Centre of Major Axis
  XCent = mean([Major_X2 , Major_X1]) ;
  YCent = mean([Major_Y2 , Major_Y1]) ;

  % Minor line equation
  Minor_Function = @(x) (Minor_Slope*(x-XCent)+YCent) ;

  % Distance from Minor Function
  DMinLine = abs(Minor_Function(X) - Y) ;

  % Finding the points that have the smallest distance from the minor equation
  Sorted_DMinLine = sort(DMinLine) ;

  % The list of points that are less than 1.~ pixel apart on the minor line
  Sorted_DMinLine(Sorted_DMinLine>1.1) = [] ; %%% It may not be formed in small forms

  % Finding points on both sides of the minor line that have the shortest distance to the center
  % Aggregation of all indexes with distance conditions of less than 1 pixel
  All_index = [] ;
  for i = 1:length(Sorted_DMinLine)
    XY_index = find(DMinLine == Sorted_DMinLine(i)) ;
    All_index = [All_index ; XY_index] ;
  end

  % Separate points on both sides of the minor
  Diff = X(All_index)-XCent ;

  % Both POS and NEG sides
  Neg = max(Diff(Diff < 0)) ;
  Pos = min(Diff(Diff > 0)) ;

  % Coordinates Minor
  Combo1 = find(Diff == Neg) ;
  Minor_X1 = X(All_index(Combo1(1))) ;
  Minor_Y1 = Y(All_index(Combo1(1))) ;
  Combo2 = find(Diff == Pos) ;
  Minor_X2 = X(All_index(Combo2(1))) ;
  Minor_Y2 = Y(All_index(Combo2(1))) ;

  % Major and Minor length
  Major = sqrt((Major_X1-Major_X2).^2 + (Major_Y1-Major_Y2).^2) ;
  Minor = sqrt((Minor_X1-Minor_X2).^2 + (Minor_Y1-Minor_Y2).^2) ;

  AspectRatioFactor = Minor/Major ;

  function [r,c] = MajorAxisinShape(X,Y,BW,D,SortD,QF)
    count = 0 ; Quality = 0 ;
    while Quality < QF
      [R,C] = find(D == SortD(end-count)) ;
      count = count + 1 ;
      % Check all indices
      for i = 1:length(R)
        % Creating a filter from indices
        Filter = poly2mask([X(R(i)) X(C(i)) X(R(i))] , [Y(R(i)) Y(C(i)) Y(R(i))],size(BW,1),size(BW,2)) ;
        % Does the created linear filter match the shape?
        Quality = nnz(logical(Filter.*BW)) / nnz(Filter) ;
        r = R(i) ; c = C(i) ;

        if Quality >= QF % If it has the necessary quality, it has succeeded in finding the major axis
          break
        end

      end
    end
  end

  end
