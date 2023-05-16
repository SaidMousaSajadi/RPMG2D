function [XMCC , YMCC , RMCC] = MCC_Finder(X,Y,BW,N,M,DMCC,Method,verbose)
  switch lower(Method)
    case 'clustering'
      % N is Threshold Point to Clustring Need 3 points to fit
      % M is Number of neighborhood
      while true
        if verbose
          disp(N)
        end
        [K , ~] = kmeans([X,Y],N) ;
##        C = {'*r' , '*b' , '*m'} ;
##        for i = 1:N
##          hold on ; plot(X(K==i),Y(K==i),C{i}) ;
##        end
        clear MCCPoints MCPoints Lside Rside
        MCCPoints = [] ;
        for i = 1:N
          ind = find(DMCC == max(DMCC(K == i))) ;
          MCPoints(i,:) = [X(ind(1)) Y(ind(1))] ;
          for j = 1:M
            Lside(j,:) = [X(ind(1)-j) Y(ind(1)-j)] ;
          end
          for j = 1:M
            Rside(j,:) = [X(ind(1)+j) Y(ind(1)+j)] ;
          end
          MCCPoints = [MCCPoints ; MCPoints ; Lside ; Rside] ;
        end
##        hold on ; plot(MCCPoints(:,1),MCCPoints(:,2),'*y')

        CentMCC = CircleFitByLandau(MCCPoints) ;
        XMCC = round(CentMCC(1)) ;
        YMCC = round(CentMCC(2)) ;
        RMCC = CentMCC(3) ;
        %     plot(XMCC,YMCC,'*y')

        LMCC = CheckCenterInGrain(BW,XMCC,YMCC) ;
        LCMCC = CheckGrainInCircle(X,Y,XMCC,YMCC,RMCC,0.90) ;

        if (LMCC == true & LCMCC == true) || N == 2 %
          break
        else
          N = N-1 ;
        end %% If

      end

    case 'peak'
      MaxLoc = LocalMaxFinder(DMCC) ;
      P = DMCC(MaxLoc) ;
##      figure ;
##      plot(DMCC) ; text(MaxLoc,P,'max') ;

      % N is Threshold Point to Clustring Need 3 points to fit
      N = length(P) ;
      Maxi = sort(P)(end:-1:1) ;

      % M is Number of neighborhood
      while true
        clear MCCPoints MCPoints Lside Rside
        MCCPoints = [] ;
        if verbose
          disp(N)
        end
        for i = 1:N
          ind = find(DMCC == Maxi(i)) ;
          MCPoints(i,:) = [X(ind(1)) Y(ind(1))] ;
          for j = 1:M
            Lside(j,:) = [X(ind(1)-j) Y(ind(1)-j)] ;
          end
          for j = 1:M
            Rside(j,:) = [X(ind(1)+j) Y(ind(1)+j)] ;
          end
          MCCPoints = [MCCPoints ; MCPoints ; Lside ; Rside] ;
        end
##        hold on ; plot(MCCPoints(:,1),MCCPoints(:,2),'*y')

        CentMCC = CircleFitByLandau(MCCPoints) ;
        XMCC = round(CentMCC(1)) ;
        YMCC = round(CentMCC(2)) ;
        RMCC = CentMCC(3) ;
        %     plot(XMCC,YMCC,'*y')

        LMCC = CheckCenterInGrain(BW,XMCC,YMCC) ;
        LCMCC = CheckGrainInCircle(X,Y,XMCC,YMCC,RMCC,0.90) ;

        if (LMCC == true && LCMCC == true) || N == 2
          break
        else
          N = N-1 ;
          Maxi(end) = [] ;
        end %% If
      end

  end %% Switch
  function MaxLoc = LocalMaxFinder(DMCC)
    Sign = sign(diff(DMCC)./1) ;
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
    MaxLoc = find((sign(DSign)<0)==1) ;
  end

end %% Function
