  function [X1 , Y1 , X2 , Y2 , X3 , Y3 , X4 , Y4 , Major , Minor , ElongationFactor] = Elongation(BW)
  % find boundary
  [Bound,Lab,~] = bwboundaries(BW) ; 
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  
  % find centre of shape
  S = regionprops(BW,'Centroid','BoundingBox','Orientation') ; % Center must be inside the shape
  
  x = [X' ; Y'] ;
  
  k = convhull(x(1,:),x(2,:));
  CH = x(:,k);
  
  E = diff(CH,1,2);           % CH edges
  T = atan2(E(2,:),E(1,:));   % angle of CH edges (used for rotation)
  T = unique(mod(T,pi/2));    % reduced to the unique set of first quadrant angles
  
  R = cos(reshape(repmat(T,2,2),2*length(T),2) + repmat([0 -pi ; pi 0]/2,length(T),1));
  
  RCH = R*CH;
  
  bsize = max(RCH,[],2) - min(RCH,[],2);
  area  = prod(reshape(bsize,2,length(bsize)/2));
  
  [a,i] = min(area);
  
  Rf    = R(2*i+[-1 0],:);   % rotated frame
  bound = Rf * CH;           % project CH on the rotated frame
  bmin  = min(bound,[],2);
  bmax  = max(bound,[],2);
  
  % compute the corner of the bounding box
  Rf = Rf';
  bb(:,1) = bmin(1)*Rf(:,1) + bmin(2)*Rf(:,2);
  bb(:,2) = bmin(1)*Rf(:,1) + bmax(2)*Rf(:,2);
  bb(:,3) = bmax(1)*Rf(:,1) + bmax(2)*Rf(:,2);
  bb(:,4) = bmax(1)*Rf(:,1) + bmin(2)*Rf(:,2);
  
  X1 = bb(1,1) ;
  Y1 = bb(2,1) ;
  X2 = bb(1,2) ;
  Y2 = bb(2,2) ;
  X3 = bb(1,3) ;
  Y3 = bb(2,3) ;
  X4 = bb(1,4) ;
  Y4 = bb(2,4) ;
  
  Axis1 = sqrt((X1-X2).^2  + (Y1-Y2).^2) ;
  Axis2 = sqrt((X2-X3).^2  + (Y2-Y3).^2) ;
  Major(Axis1 >= Axis2) = Axis1 ; Major(Axis1 < Axis2) = Axis2 ;
  Minor(Axis1 >= Axis2) = Axis2 ; Minor(Axis1 < Axis2) = Axis1 ;
  
  ElongationFactor = Minor/Major ;
  
end