function [xcent , ycent , r , XMIC , YMIC , RMIC , RoundnessFactor] = Roundness(BW,Method,TFactor,MICMethod,verbose)
  Minimum_Acceptable_Radius = 1 ;
  NeighborPixel = [15 1] ;
  % find boundary
  [Bound,~,~] = bwboundaries(BW) ;
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;

  out = cornerMyMethod(BW) ;
  XCON = out(:,1) ;
  YCON = out(:,2) ;
##  hold on ; plot(XCON,YCON,'*b')

  switch lower(Method)
    case 'withclust'
      K = evalcluster(out,2) ; % mean distance of members to center of cluster is less than 3 pixels
      index = kmeans(out,K) ;
##      for i = 1:K
##        hold on ; plot(XCON(index==i),YCON(index==i),'marker','*','Color',rand(1,3))
##      end
      r = zeros(1,K) ;
      for i = 1:K
        NeighborP = NeighborPixel(1) ;
        Flag = true ;
        while (Flag && NeighborP<=NeighborPixel(1) && NeighborP>=NeighborPixel(2))
          C{i,1} = neighborsfinder(X,Y,XCON(index==i),YCON(index==i),NeighborP) ;
          CircleProp = CircleFitByLandau([X(C{i,1}),Y(C{i,1})]) ;
          LCent(i) = CheckCenterInGrain(BW,CircleProp(1),CircleProp(2)) ;
          [LCirc(i) , Factor] = CheckCircleInGrain(X,Y,CircleProp(1),CircleProp(2),CircleProp(3),TFactor) ;
          if verbose
            clc
            disp([num2str((i/K)*100) '%'])
            disp(NeighborP)
          end
          if LCent(i) == 1 && LCirc(i) == 1
            ffactor(i) = Factor ;
            xcent(i) = CircleProp(1) ;
            ycent(i) = CircleProp(2) ;
            r(i) = CircleProp(3) ;
            Flag = false ;
          elseif LCent(i) == 1
            FFactor(NeighborP) = Factor ;
            Xcent(NeighborP) = CircleProp(1) ;
            Ycent(NeighborP) = CircleProp(2) ;
            R(NeighborP) = CircleProp(3) ;
          else
            FFactor(NeighborP) = 0 ;
            Xcent(NeighborP) = 0 ;
            Ycent(NeighborP) = 0 ;
            R(NeighborP) = 0 ;
          end
          NeighborP = NeighborP - 1 ;
        end
        if r(i) == 0
          ID = find(max(FFactor)==FFactor) ;
          ffactor(i) = FFactor(ID(1)) ;
          xcent(i) = Xcent(ID(1)) ;
          ycent(i) = Ycent(ID(1)) ;
          r(i) = R(ID(1)) ;
        end
      end
      MinPixel = Minimum_Acceptable_Radius ; % Minimum acceptable radius
      xcent(r<=MinPixel) = [] ;
      ycent(r<=MinPixel) = [] ;
      ffactor(r<=MinPixel) = [] ;
      r(r<=MinPixel) = [] ;

      MinFactor = TFactor*0.97 ; % Minimum acceptable Factor
      xcent(ffactor<=MinFactor) = [] ;
      ycent(ffactor<=MinFactor) = [] ;
      r(ffactor<=MinFactor) = [] ;
      ffactor(ffactor<=MinFactor) = [] ;

    case 'withoutclust'

      K = length(XCON) ;
      r = zeros(1,K) ;

      for i = 1:K
        NeighborP = NeighborPixel(1) ;
        Flag = true ;
        while (Flag && NeighborP<=NeighborPixel(1) && NeighborP>=NeighborPixel(2))
          C{i,1} = neighborfinder(X,Y,XCON(i),YCON(i),NeighborP) ;
          CircleProp = CircleFitByLandau([X(C{i,1}),Y(C{i,1})]) ;

          LCent(i) = CheckCenterInGrain(BW,CircleProp(1),CircleProp(2)) ;
          [LCirc(i) , Factor] = CheckCircleInGrain(X,Y,CircleProp(1),CircleProp(2),CircleProp(3),TFactor) ;

          if verbose
            clc
            disp([num2str((i/K)*100) '%'])
            disp(NeighborP)
          end

          if LCent(i) == 1 && LCirc(i) == 1
            ffactor(i) = Factor ;
            xcent(i) = CircleProp(1) ;
            ycent(i) = CircleProp(2) ;
            r(i) = CircleProp(3) ;
            Flag = false ;
          elseif LCent(i) == 1
            FFactor(NeighborP) = Factor ;
            Xcent(NeighborP) = CircleProp(1) ;
            Ycent(NeighborP) = CircleProp(2) ;
            R(NeighborP) = CircleProp(3) ;
          else
            FFactor(NeighborP) = 0 ;
            Xcent(NeighborP) = 0 ;
            Ycent(NeighborP) = 0 ;
            R(NeighborP) = 0 ;
          end
          NeighborP = NeighborP - 1 ;
        end
        if r(i) == 0
          ID = find(max(FFactor)==FFactor) ;
          ffactor(i) = FFactor(ID(1)) ;
          xcent(i) = Xcent(ID(1)) ;
          ycent(i) = Ycent(ID(1)) ;
          r(i) = R(ID(1)) ;
        end
      end

      MinPixel = Minimum_Acceptable_Radius ; % Minimum acceptable radius
      xcent(r<=MinPixel) = [] ;
      ycent(r<=MinPixel) = [] ;
      ffactor(r<=MinPixel) = [] ;
      r(r<=MinPixel) = [] ;

      MinFactor = TFactor*0.97 ; % Minimum acceptable Factor
      xcent(ffactor<=MinFactor) = [] ;
      ycent(ffactor<=MinFactor) = [] ;
      r(ffactor<=MinFactor) = [] ;
      ffactor(ffactor<=MinFactor) = [] ;

  end % switch
  S = regionprops(BW,'Centroid') ;
  Th_Clus = 5 ;
  Neighbor_Pix = 1 ;
  Th_Eff = 0.90 ;
  [XMIC , YMIC , RMIC , ~] = MIC_Finder(X,Y,S,BW,Th_Clus,Neighbor_Pix,Th_Eff,MICMethod,verbose) ;
  %% Result
  RoundnessFactor = sum(r)/(length(r)*RMIC) ;
end % function
