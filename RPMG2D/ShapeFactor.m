function [X , Y , XC , YC , Major_X1 , Major_X2 , Major_Y1 , Major_Y2 , L , ShapeFactor_Circl , ShapeFactor_Round] = ShapeFactor(BW,varargin)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Centroid','Perimeter','Area') ; % Center must be inside the shape
  
  XC = S.Centroid(1) ; 
  YC = S.Centroid(2) ; 
  
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
    
    L = sqrt((Major_X1-Major_X2).^2 + (Major_Y1-Major_Y2).^2) ;
    
    
    ShapeFactor_Circl = (4*pi*S.Area) / (S.Perimeter^2) ;
    ShapeFactor_Round = (4*S.Area) / (pi*(L^2)) ;
    
        
        
        
        
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