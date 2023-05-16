function [XMIC , YMIC , RMIC , EMIC] = MIC_Finder(X,Y,S,BW,n,m,ThresholdEff,Method,verbose)
  switch lower(Method)
    case 'centroid'
      % distance between center and boundary
      for i = 1:length(X)
        D(i,1) = pdist([S.Centroid(1) S.Centroid(2) ; X(i) Y(i)]) ;
      end

      % rearrange distance and boundary
      base_ind = find(D == max(D)) ;
      XX = [X(base_ind(1)+1:end) ; X(1:base_ind(1))] ;
      YY = [Y(base_ind(1)+1:end) ; Y(1:base_ind(1))] ;
      DD = [D(base_ind(1)+1:end) ; D(1:base_ind(1))] ;
      X = XX ; Y = YY ; D = DD ;

      %% Smooth distance (Optional)
      xD = (1:length(D))' ;
      [yD , ~ , ~ , ~] = csaps(xD,D, 0.98, xD, w=[]) ;
      D = yD ;

      %% Main Algorithm
      % Initialize N , M , TFactor
      N = n ; % Maximum number of clusters
      M = m ; % Neigber
      Tfactor = ThresholdEff ; % Quality factor

      while N>1 % for clusters: 5 4 3 2

        % Initialize x , y , r for each iteration
        x = X ; y = Y ;
        r = 0 ;

        % Clustering data with N cluster
        [K,~] = kmeans([x,y],N) ;

        % Try to estimate 15 times
        while r < 15

          % The series of points with M neighborhood and N cluster is used to fit circle
          clear MICPoints MIPoints Lside Rside
          MICPoints = [] ;
          for i = 1:N
            indx = find(D == min(D(K == i))) ;
            MIPoints(i,:) = [x(indx(1)) y(indx(1))] ;
            for j = 1:M
              Lside(j,:) = [x(indx(1)-j) y(indx(1)-j)] ;
            end
            for j = 1:M
              Rside(j,:) = [x(indx(1)+j) y(indx(1)+j)] ;
            end
            % Make series points in MICPoints
            MICPoints = [MICPoints ; MIPoints ; Lside ; Rside] ;
          end

          % Fit circle with boundary points
          CentMIC = CircleFitByLandau(MICPoints) ;
          r = r + 1 ;

          % Estimated center and radius
          rmic(r) = CentMIC(3) ;
          xmic(r) = CentMIC(1) ;
          ymic(r) = CentMIC(2) ;

          % Algorithm advance
          if verbose
            clc
            disp(["Cluster:" , num2str(N) , ", Time:" ,num2str(r)])
          end

          % Main Constraction (If the estimated center and circle are placed in the shape, save it)
          if ~isnan(rmic(r))
            LCentMIC(r) = CheckCenterInGrain(BW,xmic(r),ymic(r)) ;
            if LCentMIC(r) == true
              [LCircMIC(r) , emic(r)] = CheckCircleInGrain(x,y,xmic(r),ymic(r),rmic(r),Tfactor) ;
            else
              LCentMIC(r) = false ;
              LCircMIC(r) = false ;
              emic(r) = 0 ;
            end
          end

        end % introier while

        % Checking the results of 15 estimations
        if ~isempty(rmic)
          ne = 2*MyNormalize(emic,"range") + MyNormalize(rmic,"range") ;
          mic_index = find(ne == max(ne)) ; % best
          rMIC(N) = rmic(mic_index(1)) ;
          xMIC(N) = xmic(mic_index(1)) ;
          yMIC(N) = ymic(mic_index(1)) ;
          eMIC(N) = emic(mic_index(1)) ;
        else
          rMIC(N) = nan ;
          xMIC(N) = nan ;
          yMIC(N) = nan ;
          eMIC(N) = nan ;
        end

        % reduce of the clusters and check the algorithm again
        N = N-1 ;

      end % Main Loop
      % Get XMIC,YMIC,RMIC,EMIC
      Ne = 2*MyNormalize(eMIC,"range") + MyNormalize(rMIC,'range') ;
      MIC_index = find(Ne == max(Ne)) ; % best
      RMIC = rMIC(MIC_index(1)) ;
      EMIC = eMIC(MIC_index(1)) ;
      XMIC = xMIC(MIC_index(1)) ;
      YMIC = yMIC(MIC_index(1)) ;
    case 'skeleton'

      % Make Skeleton
      bwth = bwmorph(BW,'thin',Inf) ;
      filter = [1 1 1 ;
      1 0 1 ;
      1 1 1] ; % filter
      bwdisc1 = bwth & ~(bwth & conv2(double(bwth),filter,'same') > 1) ;
      ##figure ; imshow(bwdisc1)
      bwdisc2 = bwth & ~(bwth & conv2(double(bwth),filter,'same') > 2) ;
      ##figure ; imshow(bwdisc2)

      bwcc1 = bwconncomp(bwdisc1) ;
      bwcc2 = bwconncomp(bwdisc2) ;

      XC = [] ;
      YC = [] ;

      if bwcc2.NumObjects > bwcc1.NumObjects
        for j = 1:bwcc2.NumObjects
          A = cell2mat(bwcc1.PixelIdxList(1,:)) ;
          if isempty(intersect(bwcc2.PixelIdxList{1,j} , A))
            [YCr,XCr] = ind2sub(size(bwdisc2),bwcc2.PixelIdxList{1,j}) ;
            YC = [YC ; YCr] ; XC = [XC ; XCr] ;
          end
        end
      else
        for j = 1:bwcc2.NumObjects
          [YCr,XCr] = ind2sub(size(bwdisc2),bwcc2.PixelIdxList{1,j}) ;
          YC = [YC ; YCr] ; XC = [XC ; XCr] ;
        end
      end
      % hold on ; plot(XC,YC,'or')
      OrgX = X ; OrgY = Y ;
      for k = 1:length(XC)
        X = OrgX ; Y = OrgY ;
        SXC = XC(k) ;
        SYC = YC(k) ;
        for i = 1:length(X)
          D(i,1) = pdist([SXC SYC ; X(i) Y(i)]) ;
        end

        % rearrange distance and boundary
        base_ind = find(D == max(D)) ;
        XX = [X(base_ind(1)+1:end) ; X(1:base_ind(1))] ;
        YY = [Y(base_ind(1)+1:end) ; Y(1:base_ind(1))] ;
        DD = [D(base_ind(1)+1:end) ; D(1:base_ind(1))] ;
        X = XX ; Y = YY ; D = DD ;

        %% Smooth distance (Optional)
        xD = (1:length(D))' ;
        [yD , ~ , ~ , ~] = csaps(xD,D, 0.98, xD, w=[]) ;
        D = yD ;

        %% Main Algorithm
        % Initialize N , M , TFactor
        N = n ; % Maximum number of clusters
        M = m ; % Neigber
        Tfactor = ThresholdEff ; % Quality factor

        while N>1 % for clusters: 5 4 3 2

          % Initialize x , y , r for each iteration
          x = X ; y = Y ;
          r = 0 ;

          % Clustering data with N cluster
          [K,~] = kmeans([x,y],N) ;

          % Try to estimate 15 times
          while r < 15

            % The series of points with M neighborhood and N cluster is used to fit circle
            clear MICPoints MIPoints Lside Rside
            MICPoints = [] ;
            for i = 1:N
              indx = find(D == min(D(K == i))) ;
              MIPoints(i,:) = [x(indx(1)) y(indx(1))] ;
              for j = 1:M
                Lside(j,:) = [x(indx(1)-j) y(indx(1)-j)] ;
              end
              for j = 1:M
                Rside(j,:) = [x(indx(1)+j) y(indx(1)+j)] ;
              end
              % Make series points in MICPoints
              MICPoints = [MICPoints ; MIPoints ; Lside ; Rside] ;
            end

            % Fit circle with boundary points
            CentMIC = CircleFitByLandau(MICPoints) ;
            r = r + 1 ;

            % Estimated center and radius
            rmic(r) = CentMIC(3) ;
            xmic(r) = CentMIC(1) ;
            ymic(r) = CentMIC(2) ;

            % Algorithm advance
            if verbose
              clc
              disp(["Cluster:" , num2str(N) , ", Time:" ,num2str(r)])
            end

            % Main Constraction (If the estimated center and circle are placed in the shape, save it)
            if ~isnan(rmic(r))
              LCentMIC(r) = CheckCenterInGrain(BW,xmic(r),ymic(r)) ;
              if LCentMIC(r) == true
                [LCircMIC(r) , emic(r)] = CheckCircleInGrain(x,y,xmic(r),ymic(r),rmic(r),Tfactor) ;
              else
                LCentMIC(r) = false ;
                LCircMIC(r) = false ;
                emic(r) = 0 ;
              end
            end

          end % introier while

          % Checking the results of 15 estimations
          if ~isempty(rmic)
            ne = 2*MyNormalize(emic,"range") + MyNormalize(rmic,"range") ;
            if ~all(isnan(ne))
              mic_index = find(ne == max(ne)) ; % best
              rMIC(N) = rmic(mic_index(1)) ;
              xMIC(N) = xmic(mic_index(1)) ;
              yMIC(N) = ymic(mic_index(1)) ;
              eMIC(N) = emic(mic_index(1)) ;
            else
              rMIC(N) = nan ;
              xMIC(N) = nan ;
              yMIC(N) = nan ;
              eMIC(N) = nan ;
            end
          else
            rMIC(N) = nan ;
            xMIC(N) = nan ;
            yMIC(N) = nan ;
            eMIC(N) = nan ;
          end

          % reduce of the clusters and check the algorithm again
          N = N-1 ;

        end % Main Loop
        if ~isempty(rMIC)
          Ne_Skel = 2*MyNormalize(eMIC,"range") + MyNormalize(rMIC,'range') ;
          if ~all(isnan(Ne_Skel))
            MIC_index_Skel = find(Ne_Skel == max(Ne_Skel)) ; % best
            RMIC_Skel = rMIC(MIC_index_Skel(1)) ;
            EMIC_Skel = eMIC(MIC_index_Skel(1)) ;
            XMIC_Skel = xMIC(MIC_index_Skel(1)) ;
            YMIC_Skel = yMIC(MIC_index_Skel(1)) ;
            RMIC_Skel_Cal(k) = RMIC_Skel ;
            EMIC_Skel_Cal(k) = EMIC_Skel ;
            XMIC_Skel_Cal(k) = XMIC_Skel ;
            YMIC_Skel_Cal(k) = YMIC_Skel ;
          else
            RMIC_Skel_Cal(k) = nan ;
            EMIC_Skel_Cal(k) = nan ;
            XMIC_Skel_Cal(k) = nan ;
            YMIC_Skel_Cal(k) = nan ;
          end
        else
          RMIC_Skel_Cal(k) = nan ;
          EMIC_Skel_Cal(k) = nan ;
          XMIC_Skel_Cal(k) = nan ;
          YMIC_Skel_Cal(k) = nan ;
        end
      end
      % Get XMIC,YMIC,RMIC,EMIC
      if ~isempty(RMIC_Skel_Cal)
        Ne_Skel_Cal = 2*MyNormalize(EMIC_Skel_Cal,"range") + MyNormalize(RMIC_Skel_Cal,'range') ;
        if ~all(isnan(Ne_Skel_Cal))
          MIC_index_Skel_Cal = find(Ne_Skel_Cal == max(Ne_Skel_Cal)) ; % best
          RMIC = RMIC_Skel_Cal(MIC_index_Skel_Cal(1)) ;
          EMIC = EMIC_Skel_Cal(MIC_index_Skel_Cal(1)) ;
          XMIC = XMIC_Skel_Cal(MIC_index_Skel_Cal(1)) ;
          YMIC = YMIC_Skel_Cal(MIC_index_Skel_Cal(1)) ;
        else
          RMIC = nan ;
          EMIC = nan ;
          XMIC = nan ;
          YMIC = nan ;
        end
      else
        RMIC = nan ;
        EMIC = nan ;
        XMIC = nan ;
        YMIC = nan ;
      end
  endswitch
  end
