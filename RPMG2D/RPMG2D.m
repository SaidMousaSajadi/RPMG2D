close all ; clc ; clear all ;


graphics_toolkit qt % or fltk(PLC/old) , qt(apps) , gnuplot(sites)
pkg load image % images
pkg load statistics % pdist
pkg load signal % findpeaks
pkg load splines % smoothing lines
pkg load ltfat % normalize
pkg load optim
pkg load io

global NewD
global K
root.theta = 0:0.01:2*pi ;
root.FlagSteps = false ;
root.Shape = 0 ;
root.FlagGeneShape = false ;
root.FlagGenePM = false ;

%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
function Theta = AngleFinder(x0,y0,X,Y)
  Theta = atan2(y0 - Y, x0 - X) * (180 / pi) ;
end
function [XC,YC,W,X,Y,D,T] = InitialValue(BW)
  [Bound,~,~] = bwboundaries(BW) ;
  X = Bound{1,1}(:,2) ;
  Y = Bound{1,1}(:,1) ;
  S = regionprops(BW,'Centroid') ;
  XC = S.Centroid(1) ;
  YC = S.Centroid(2) ;
  T = AngleFinder(XC,YC,X,Y) ;
  for i = 1:length(X)
    D(i) = pdist([XC YC ; X(i) Y(i)]) ;
  end
  D = D' ;
  t = AngleFinder(S.Centroid(1),S.Centroid(2),X,Y) ;
  W = 1 ;
end
function [x,y,d,t] = ReArrange(X,Y,D,T)
  % rearrange distance and boundary
  base_ind = find(D == max(D)) ;
  x = [X(base_ind(1)+1:end) ; X(1:base_ind(1))]  ;
  y = [Y(base_ind(1)+1:end) ; Y(1:base_ind(1))]  ;
  d = [D(base_ind(1)+1:end) ; D(1:base_ind(1))] ;
  t = [T(base_ind(1)+1:end) ; T(1:base_ind(1))] ;
end
function D = Smoother(D,Factor)
  xD = (1:length(D))' ;
  [D , ~ , ~ , ~] = csaps(xD,D,Factor, xD, w=[]) ;
end
function [X,Y,XC,YC] = BoundGenerator(x0,y0,D,W,t,Rotate)
  % Rotate degree
  X=x0+((W * D).*cosd(t+180+Rotate));
  Y=y0+((W * D).*sind(t+180+Rotate));
  % Connect Start Point to End Point
  m = (Y(1)-Y(end))./(X(1)-X(end)) ;
  f = @(x) m*(x-X(1)) + Y(1) ;
  xD = linspace(X(1),X(end),20) ;
  yD = f(xD) ;
  X = [X ; xD'] ;
  Y = [Y ; yD'] ;
  XC = x0 ;
  YC = y0 ;
  if any(X<=1) || any(Y<=1) % replace shape
    X = X+abs(min(X))+5 ;
    Y = Y+abs(min(Y))+5 ;
    XC = x0+abs(min(X))+5 ;
    YC = y0+abs(min(Y))+5 ;
  end
end
function [Lower , Upper , Guess] = BoundFinder(NewD,K,index)
  Lower = [] ;
  Upper = [] ;
  % Coefficient Bound
  for i = 1:K
    Lower = [Lower , max((NewD(index==i)))*0.1 ] ;
    Upper = [Upper , max((NewD(index==i)))*1.5 ] ;
  end
  % Mu Bound
  for i = 1:K
    mu = find(max(NewD(index==i))==NewD) ;
    Lower = [Lower , mu-10 ] ;
    Upper = [Upper , mu+10] ;
  end
  % Sigma Bound
  for i = 1:K
    Lower = [Lower , 0 ] ;
    Upper = [Upper , 100] ;
  end
  Guess = mean([Lower ; Upper]) ;
end
function [NewX,NewY,CC] = XY2BWMatrix(ScaledX,ScaledY,Size)
  BOW = zeros(Size) ;
  BOW(sub2ind(Size,round(ScaledY),round(ScaledX))) = 1 ;
  CC = imdilate(BOW,[0 1 0 ; 1 1 1 ; 0 1 0]) ;
  CC = imclose(CC,[0 1 0 ; 1 1 1 ; 0 1 0]) ;
  CC = bwmorph(CC,'skel',Inf);
  CC = imfill(CC,'holes');
  [NewBound,~,~] = bwboundaries(CC) ;
  NewX = NewBound{1,1}(:,2) ;
  NewY = NewBound{1,1}(:,1) ;
end
function [AINumber_S1 , AINumber_S2 , AINumber_S3] = ModelSelector(AINumber,NewAINumber_index,Method)
  switch lower(Method)
    case 'vcm'
      % Variant Coefficient model(VCM) - moderate similar
      AINumber_S1 = [AINumber(NewAINumber_index(1,:),1) AINumber(:,2) AINumber(:,3)] ;
      AINumber_S2 = [AINumber(NewAINumber_index(2,:),1) AINumber(:,2) AINumber(:,3)] ;
      AINumber_S3 = [AINumber(NewAINumber_index(3,:),1) AINumber(:,2) AINumber(:,3)] ;
    case 'vmm'
      % Variant Mu model(VMM) - poorly similar
      AINumber_S1 = [AINumber(:,1) AINumber(NewAINumber_index(1,:),2) AINumber(:,3)] ;
      AINumber_S2 = [AINumber(:,1) AINumber(NewAINumber_index(2,:),2) AINumber(:,3)] ;
      AINumber_S3 = [AINumber(:,1) AINumber(NewAINumber_index(3,:),2) AINumber(:,3)] ;
    case 'vsm'
      % Variant Sigma model(VSM) - very similar
      AINumber_S1 = [AINumber(:,1) AINumber(:,2) AINumber(NewAINumber_index(1,:),3)] ;
      AINumber_S2 = [AINumber(:,1) AINumber(:,2) AINumber(NewAINumber_index(2,:),3)] ;
      AINumber_S3 = [AINumber(:,1) AINumber(:,2) AINumber(NewAINumber_index(3,:),3)] ;
  endswitch
  AINumber_S1 = (5*AINumber+AINumber_S1)./6 ;
  AINumber_S2 = (5*AINumber+AINumber_S2)./6 ;
  AINumber_S3 = (5*AINumber+AINumber_S3)./6 ;
end
function NewAINumber_index = FindIndex2NewAINumbers(K)
  if K == 1 % Circle
    NewAINumber_index = ones(3,1) ;
  elseif K == 2  % Line
    NewAINumber_index = [1 , 2 ; 1 , 2 ; 1 , 2] ;
  elseif K == 3 % Crossed lines
    NewAINumber_index = [1 , 2 , 3 ; 1 , 2 , 3 ; 1 , 2 , 3] ;
  elseif K == 4 % Triangle
    NewAINumber_index = [1 , 2 , 3 , 4 ; 1 , 3 , 2 , 4 ; 1 , 2 , 3 , 4] ;
  else
    NewAINumber_index = [ones(3,1) perms_k([1+1:K-1],3) K*ones(3,1)] ;
  end
end
function TableHandle = Table_Update(TableHandle,ShapeHandle,Value,row)
  LocalMatrix = get(TableHandle,'Data') ;
  LocalMatrix(row,ShapeHandle+1) = Value ;
  set(TableHandle,'Data',LocalMatrix) ;
end
function SturctToJSON(Cell,NameFile)
  Shapes_JSON = jsonencode(Cell,'PrettyPrint',true) ;
  oid = fopen(NameFile,'wt');
  fprintf(oid, Shapes_JSON);
  fclose(oid);
end
function ShowAngularity(Ax,BW,xc,yc,r,theta,xfar,yfar,xcirc,ycirc)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  xcircle = xc + r*cos(theta) ;
  ycircle = yc + r*sin(theta) ;
  plot(Ax,xcircle,ycircle,"r")
  plot(Ax,[xc xfar],[yc yfar],'-b','LineWidth',1.25)
  plot(Ax,[xc xcirc],[yc ycirc],'-b','LineWidth',1.25)
  plot(Ax,[xfar xcirc],[yfar ycirc],'-b','LineWidth',1.25)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowAspectRatio(Ax,BW,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,[majorx1 , majorx2],[majory1 , majory2],'-b','LineWidth',1.25) ;
  plot(Ax,[minorx1 , minorx2],[minory1 , minory2],'-r','LineWidth',1.25) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowCircularity(Ax,BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,theta)
  xmic_circ = xmic + rmic*cos(theta) ;
  ymic_circ = ymic + rmic*sin(theta) ;
  xmcc_circ = xmcc + rmcc*cos(theta) ;
  ymcc_circ = ymcc + rmcc*sin(theta) ;
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,xmic,ymic,'*b') ;
  plot(Ax,xmic_circ,ymic_circ,'b','LineWidth',1.75) ;
  plot(Ax,xmcc,ymcc,'*r') ;
  plot(Ax,xmcc_circ,ymcc_circ,'r','LineWidth',1.75) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowMZC(Ax,BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,theta)
  xmcc_circ = xmcc + rmcc*cos(theta) ;
  ymcc_circ = ymcc + rmcc*sin(theta) ;
  xmic_circ = xmic + rmic*cos(theta) ;
  ymic_circ = ymic + rmic*sin(theta) ;
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  patch(Ax,[xmic_circ , xmcc_circ] , [ymic_circ , ymcc_circ] , 'g','linestyle','none','FaceAlpha' , 0.3)
  plot(Ax,xmic_circ,ymic_circ,'-b','LineWidth',1.5)
  plot(Ax,xmcc_circ,ymcc_circ,'-b','LineWidth',1.5)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end

function ShowLSC(Ax,BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,theta)
  rlsc = mean([rmic,rmcc]) ;
  xlsc = mean([xmic,xmcc]) ;
  ylsc = mean([ymic,ymcc]) ;
  xlsc_circ = xlsc + rlsc*cos(theta) ;
  ylsc_circ = ylsc + rlsc*sin(theta) ;
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,xlsc_circ,ylsc_circ,'m','LineWidth',1.75) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowCompactness(Ax,BW,x,y)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,x,y,'b','LineWidth',1.25)
  text(Ax,mean(x),mean(y),'\Omega','Color','w','Interpreter','tex','FontSize',11) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowConvexity(Ax,BW,x,y,ConvexBoundary)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,x,y,"-b",'LineWidth',1.25)
  plot(Ax,ConvexBoundary(:,1),ConvexBoundary(:,2),"-r",'LineWidth',1.25)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowEccentricity(Ax,BW,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,Elli_Bound(1,:),Elli_Bound(2,:),'r','LineWidth',1.25) ;
  plot(Ax,[X1 X2],[Y1 Y2],'-g')
  plot(Ax,[x1 x2],[y1 y2],'-m')
  plot(Ax,[xc onMajorx] , [yc onMajory],'g','LineWidth',2.5)
  plot(Ax,[xc onMinorx] , [yc onMinory],'m','LineWidth',2.5)
  plot(Ax,[onMajorx onMinorx],[onMajory onMinory],'c','LineWidth',2.5)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowElliptic(Ax,BW,x,y,ellix,elliy)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,x,y,'b','LineWidth',1.25)
  plot(Ax,ellix,elliy,'r','LineWidth',1.25)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowElongation(Ax,BW,x1,x2,x3,x4,y1,y2,y3,y4)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,[x1,x2,x3,x4,x1],[y1,y2,y3,y4,y1],'-g','LineWidth',1.25)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowRegularity(Ax,BW,xc,yc,x,y)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,x,y,'b','LineWidth',1.25)
  plot(Ax,xc,yc,'r','LineWidth',1.25)
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowRoundness(Ax,BW,xcent,ycent,r,XMIC,YMIC,RMIC,theta)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  for i = 1:length(r)
    xcirc = xcent(i) + r(i)*cos(theta) ;
    ycirc = ycent(i) + r(i)*sin(theta) ;
    plot(Ax,xcirc,ycirc,'g','LineWidth',1.5) ;
  end
  xcircMIC = XMIC + RMIC*cos(theta) ;
  ycircMIC = YMIC + RMIC*sin(theta) ;
  plot(Ax,xcircMIC,ycircMIC,'b','LineWidth',2) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowShapeFactor(Ax,BW,x,y,xc,yc,X1,X2,Y1,Y2)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  plot(Ax,x,y,'-b','LineWidth',1.25)
  plot(Ax,[X1 X2],[Y1 Y2],'-g','LineWidth',1.25)
  text(Ax,xc,yc,'\Omega','Color','w','Interpreter','tex','FontSize',11) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowSolidity(Ax,BW,x,y,xc,yc)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  patch(Ax,xc,yc,'r','edgecolor','none')
  patch(Ax,x,y,'k','edgecolor','none')
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function ShowSphericity(Ax,BW,xc,yc,rmcc,theta)
  imshow(~BW,'parent',Ax)
  hold(Ax,'on') ;
  Xcirc = xc + rmcc*cos(theta) ;
  Ycirc = yc + rmcc*sin(theta) ;
  plot(Ax,Xcirc,Ycirc,'r','LineWidth',1.5) ;
  ID = randi(length(Xcirc)) ;
  plot(Ax,[xc Xcirc(ID)],[yc Ycirc(ID)],'r','LineWidth',1.5) ;
  axis(Ax,'tight') ;
  hold(Ax,'off') ;
end
function Processing(STR,obj)
  if STR == 'e'
    set(obj,'string',' Error Occurred !!!')
    set(obj,'backgroundcolor',[1 0.785 0.785])
  elseif STR == 'p'
    set(obj,'string',' In Processing ...')
    set(obj,'backgroundcolor',[1 0.863 0.666])
  else
    set(obj,'string',' Ready')
    set(obj,'backgroundcolor',[0.785 1 0.785])
  end
end
function [Center,R] = RandomShapeGenerator(Dim,R_min,R_max,Bound)

  %% Grids
  dx = R_min/50 ;
  x = 0:dx:Dim(1) ;
  y = 0:dx:Dim(2) ;
  [x,y] = meshgrid(x,y) ;
  Gir = [x(:) y(:)] ; % initial grids
  dR = R_max - R_min ;
  Cc = bsxfun(@times,Dim,[0 0; 1 0; 1 1; 0 1]) ; % corner vertices
  % Remove All Grids Less Than R_min From Boundary
  G_max = Gir+R_min+1E-12 ;
  G_min = Gir-R_min-1E-12 ;
  chk_in = bsxfun(@le,G_max,Dim) & bsxfun(@ge,G_min,[0 0]) ;
  chk_in = sum(chk_in,2)==2 ;
  Gir(~chk_in,:) = [] ; % initial eligible grids
  Center = [] ; R = [] ;
  Ng = size(Gir,1) ; % Number Of Grids
  flag = true ;
  n = 0 ; Count = 0; m = 0 ; % initialize
  while ~isempty(Gir) && Count<1E4
    n = n+1 ;
    % New circle
    if flag && (Count>600 || size(Gir,1)<0.95*Ng)
      flag = false ;
      Rg = R_max*ones(size(Gir,1),1) ;
    end
    i = [] ;
    % Find Center Location
    if Count<=600 && flag
      X_new = Dim.*rand(1,2) ;   % centroid
    else
      i = randi(size(Gir,1)) ;
      X_new = Gir(i,:)+(dx/2)*(2*rand(1,2)-1) ;
      X_new = min(max(X_new,[0 0]),Dim) ;
      if Count>1E3
        Rg(:) = max(0.95*Rg,R_min) ;
      end
    end
    % Find Radius
    if isempty(i)
      R_new = dR*rand(1)+R_min ; % radius
    else
      R_new = (Rg(i)-R_min)*rand(1)+R_min ;
    end
    % Check if the circle fits inside the rectangle when Bound=true
    if Bound
      X_new_max = X_new+R_new+1E-12 ;
      X_new_min = X_new-R_new-1E-12 ;
      chk_in = X_new_max<=Dim & X_new_min>=[0 0] ;
      if sum(chk_in)<2
        Count = Count+1 ;
        continue
      end
    end
    % Does the new circle intersect with any other circles?
    if ~isempty(Center)
      d_in = sqrt(sum(bsxfun(@minus,Center,X_new).^2,2)) ;
      id = d_in<(R+R_new) ;
      if sum(id)>0
        Count = Count+1 ;
        if ~isempty(i)
          Rg(i) = min(0.95*Rg(i),min(0.99*(R_new+dx/2),R_max)) ;
          Rg(i) = max(Rg(i),R_min) ;
        end
        continue
      end
    end
    % Accept new circle and Update Grids
    Count = 0 ;
    m = m+1 ;
    Center = cat(1,Center,X_new) ;
    R = cat(1,R,R_new) ;
    [Gir,id] = update_grid(X_new,R_new,Gir,R_min) ; % remaining grids
    %
    if ~flag
      Rg(id) = [] ;
    end
  end
end
function [Gir,id] = update_grid(X_new,R_new,Gir,R_min)
    % Remove grid points within R_new+R_min units of new circle
    D = sum(bsxfun(@minus,Gir,X_new).^2,2) ;
    id = D<(R_new+R_min+1E-12)^2 ;
    Gir(id,:) = [] ;
end
function [Center,R] = SampleLineSegment(Xa,Xb,Ra,Rb,R_min,R_max)
    % Place circles along line segment between points Xa and Xb
    r = Xb-Xa ;
    L = norm(r) ;
    r = r/norm(L) ;
    dR = R_max-R_min ;
    Center = Xa ; R = Ra ; % initialize
    while true
        R_new = dR*rand(1)+R_min ;
        C_new = Center(end,:)+r*(R(end)+R_new+R_max*rand(1)) ;
        % Distance
        D = L - norm(C_new + r*(R_new+Rb) - Xa) ; % will there be enough space left for the end point with radius Rb?
        if D < 2*(R_min+1E-12)
            Center = cat(1,Center,Xb) ;
            R = cat(1,R,Rb) ;
            break
        else
            Center = cat(1,Center,C_new) ;
            R = cat(1,R,R_new) ;
        end
    end
end
function R = randbtw(L,H)
    R = L+(H-L).*rand ;
end
function I = Nearst(A,Fun1,x)
    Diff = abs(Fun1-A) ;
    index = find(Diff==min(Diff)) ;
    I = x(index(1)) ;
end
function [R1,R5,R16,R50,R84,R95,R100,a,X,Fun,Category,XForFit,YForFit,PhiValue] = SortRequest(Method)
    Mean = 7 ; % Peak of Maximum Phi
    x = -40:0.01:Mean ;
    switch lower(Method)
        case {'verywell' , 1}
            a = randbtw(0.0192,0.7572) ;
        case {'well' , 2}
            a = randbtw(0.7572,1.0821) ;
        case {'moderatelywell' , 3}
            a = randbtw(1.0821,1.5108) ;
        case {'moderately' , 4}
            a = randbtw(1.5108,2.1636) ;
        case {'poorly' , 5}
            a = randbtw(2.1636,4.3303) ;
        case {'verypoorly' , 6}
            a = randbtw(4.3303,4.9820) ;
    end
    Fun1 = 100*exp(-((x-Mean)/a).^2) ;
    R1 = Nearst(1,Fun1,x) ;
    R5 = Nearst(5,Fun1,x) ;
    R16 = Nearst(16,Fun1,x) ;
    R50 = Nearst(50,Fun1,x) ;
    R84 = Nearst(84,Fun1,x) ;
    R95 = Nearst(95,Fun1,x) ;
    R100 = Nearst(100,Fun1,x) ;

    X = R1:0.01:R100 ;
    Fun = 100*exp(-((X-R100)/a).^2) ; % Cumulative Frequency Vector
    Category = sort(linspace(R100,R1-1,20)) ;
    for i = 1:length(Category)
        CFC(i) = Nearst(Category(i),X,Fun) ;
    end
    XForFit = (Category(2:end) + Category(1:end-1))/2 ;
    YForFit = diff(CFC) ; % Cumulative Frequency Category

    PhiValue = (R84-R16)/4 + (R95-R5)/6.6 ;
end
function [r] = GenerateKthRadius(xx,yy,k)
    % From the frequency function, we extract the k numbers phi and convert it to radius and check radius
    coef = sum(yy)/k ;
    newyy = yy/coef ;
    NOEC = round(newyy) ;  %Number of each category
    phi = [] ;
    for i = 1:length(NOEC)
        if NOEC(i) ~= 0
            phi = [phi xx(i)*ones(1,NOEC(i))] ;
        end
    end
    % phi = -log2(d) > phi = log2(d^-1) > 2^phi = d^-1 > d = 1/(2*phi) > 2r = 1/(2*phi) > r[mm] = 1/(2*(2^phi))
    r = 1./(2*(2.^phi)) ; %[mm]
    r = r*10 ;
    r = r*37.7952755906 ; % to pixel
    r = -sort(-r) ;
end
function [X,Y] = ShapePainter(x,y,r,theta,Rotate)
    for i = 1:length(x)
        X(i,:)=x(i)+r(i,:).*cosd(theta(i,:)+180+Rotate(i,1));
        Y(i,:)=y(i)+r(i,:).*sind(theta(i,:)+180+Rotate(i,1));
    end
end






function Exiter(obj)
  h = guidata(obj) ; % get handles
  close(h.Fig)
end
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%           %%%%%%%%%%%%%%%%%%%%%%%%

% Main
function Update_UI(obj,init = false)
  h = guidata(obj) ; % get handles
  switch (gcbo)
    case {h.Open}
      Processing('p',h.Process)
      [h.FileName, h.FilePath, h.FileIndex] = uigetfile({"*.jpeg;*.jpg;*.tiff;*.tif;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group";"*.tif ; *.tiff", "Tagged Image File Format"}) ;
      if (h.FileIndex) ~= 0
        NameSpl = strsplit(h.FileName,".") ;
        switch NameSpl{1,end}
          case {'jpg','jpeg','png','tif','tiff'}
            [h.IM , h.Map]= imread([h.FilePath h.FileName]) ;
            [h.BW , h.Map] = Image2Binary(h.IM , h.Map) ;
            h.BW = imfill(~h.BW,'holes');
            h.BW = imresize(h.BW,[200 NaN]) ; % for decrease compute process
            [XC,YC,W,X,Y,D,T] = InitialValue(h.BW) ;
            init_Org.XC = XC ; init_Org.YC = YC ; init_Org.W = W ; init_Org.X = X ; init_Org.Y = Y ; init_Org.D = D ; init_Org.Theta = T ;
            h.init_Org = init_Org ;
            imshow(~h.BW,'parent',h.Ax_Org)
            imshow(~h.BW,'parent',h.Ax_Ana)
            cla(h.Ax_Gen1)
            cla(h.Ax_Gen2)
            cla(h.Ax_Gen3)
            set(h.Table,'Data',h.TableData) ;
            h.FlagSteps = true ;
            h.FlagGeneShape = false ;
            h.FlagGenePM = false ;
            guidata(gcf,h) % update handles
        end
      end
      Processing('r',h.Process)

    case {h.SaveImgAn}
      Processing('p',h.Process)
      [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group"}) ;
      if (FileIndex) ~= 0
        F = getframe(h.Ax_Ana);
        Image = frame2im(F);
        imwrite(Image, [FilePath ,  FileName])
      end
      Processing('r',h.Process)

    case {h.SaveTable}
      Processing('p',h.Process)
      [FileName, FilePath, FileIndex] = uiputfile({"*.csv","Comma-Separated Values"}) ;
      if (FileIndex) ~= 0
        mycsvwrite([FilePath ,  FileName],get(h.Table,'RowName'),get(h.Table,'ColumnName'),get(h.Table,'Data')) ;
      end
      Processing('r',h.Process)

    case {h.SaveImgRe}
      Processing('p',h.Process)
      if h.FlagGeneShape
        [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group"}) ;
        if (FileIndex) ~= 0
          imshow(~h.NewBWShow,'parent',h.Ax_Ana)
          F = getframe(h.Ax_Ana);
          Image = frame2im(F);
          imwrite(Image, [FilePath ,  FileName])
        end
      end
      Processing('r',h.Process)

    case {h.SaveImgJSON}
      Processing('p',h.Process)
      if h.FlagGeneShape
        [FileName, FilePath, FileIndex] = uiputfile({"*.json","JavaScript Object Notation"}) ;
        if (FileIndex) ~= 0
          SturctToJSON(h.Org,[FilePath FileName])
        end
      elseif h.FlagSteps
        [FileName, FilePath, FileIndex] = uiputfile({"*.json","JavaScript Object Notation"}) ;
        if (FileIndex) ~= 0
          SturctToJSON(h.init_Org,[FilePath FileName])
        end
      else
      end
      Processing('r',h.Process)

    case {h.Radio0}
      Processing('p',h.Process)
      set (h.Radio0, "value", 1);
      set (h.Radio1, "value", 0);
      set (h.Radio2, "value", 0);
      set (h.Radio3, "value", 0);
      h.Shape = 0 ;
      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Radio1}
      Processing('p',h.Process)
      set (h.Radio0, "value", 0);
      set (h.Radio1, "value", 1);
      set (h.Radio2, "value", 0);
      set (h.Radio3, "value", 0);
      h.Shape = 1 ;
      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Radio2}
      Processing('p',h.Process)
      set (h.Radio0, "value", 0);
      set (h.Radio1, "value", 0);
      set (h.Radio2, "value", 1);
      set (h.Radio3, "value", 0);
      h.Shape = 2 ;
      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Radio3}
      Processing('p',h.Process)
      set (h.Radio0, "value", 0);
      set (h.Radio1, "value", 0);
      set (h.Radio2, "value", 0);
      set (h.Radio3, "value", 1);
      h.Shape = 3 ;
      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Angularity_Cent}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.BW,'Centroid') ;
          ShowAngularity(h.Ax_Ana,h.BW,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.NewShape1,'Centroid') ;
          ShowAngularity(h.Ax_Gen1,h.NewShape1,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.NewShape2,'Centroid') ;
          ShowAngularity(h.Ax_Gen2,h.NewShape2,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.NewShape3,'Centroid') ;
          ShowAngularity(h.Ax_Gen3,h.NewShape3,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Angularity_Skel}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.BW,'Skeleton') ;
          ShowAngularity(h.Ax_Ana,h.BW,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.NewShape1,'Skeleton') ;
          ShowAngularity(h.Ax_Gen1,h.NewShape1,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.NewShape2,'Skeleton') ;
          ShowAngularity(h.Ax_Gen2,h.NewShape2,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,r,eff,xfar,yfar,X,xcirc,ycirc,alpha,AngularFactor] = Angularity(h.NewShape3,'Skeleton') ;
          ShowAngularity(h.Ax_Gen3,h.NewShape3,xc,yc,r,h.theta,xfar,yfar,xcirc,ycirc)
          h.Table = Table_Update(h.Table,h.Shape,AngularFactor,1) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.AspectRatio_Outside}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.BW) ;
          ShowAspectRatio(h.Ax_Ana,h.BW,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.NewShape1) ;
          ShowAspectRatio(h.Ax_Gen1,h.NewShape1,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.NewShape2) ;
          ShowAspectRatio(h.Ax_Gen2,h.NewShape2,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.NewShape3) ;
          ShowAspectRatio(h.Ax_Gen3,h.NewShape3,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.AspectRatio_Inside}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.BW,'Inside') ;
          ShowAspectRatio(h.Ax_Ana,h.BW,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.NewShape1,'Inside') ;
          ShowAspectRatio(h.Ax_Gen1,h.NewShape1,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.NewShape2,'Inside') ;
          ShowAspectRatio(h.Ax_Gen2,h.NewShape2,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [majorx1,majorx2,majory1,majory2,minorx1,minory1,minorx2,minory2,major,minor,AspectFactor] = AspectRatio(h.NewShape3,'Inside') ;
          ShowAspectRatio(h.Ax_Gen3,h.NewShape3,majorx1,majorx2,majory1,majory2,minorx1,minorx2,minory1,minory2)
          h.Table = Table_Update(h.Table,h.Shape,AspectFactor,2) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Circular_Cent_Clus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.BW,'Centroid','Clustering') ;
          ShowCircularity(h.Ax_Ana,h.BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape1,'Centroid','Clustering') ;
          ShowCircularity(h.Ax_Gen1,h.NewShape1,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape2,'Centroid','Clustering') ;
          ShowCircularity(h.Ax_Gen2,h.NewShape2,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape3,'Centroid','Clustering') ;
          ShowCircularity(h.Ax_Gen3,h.NewShape3,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Circular_Cent_Peak}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.BW,'Centroid','Peak') ;
          ShowCircularity(h.Ax_Ana,h.BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape1,'Centroid','Peak') ;
          ShowCircularity(h.Ax_Gen1,h.NewShape1,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape2,'Centroid','Peak') ;
          ShowCircularity(h.Ax_Gen2,h.NewShape2,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape3,'Centroid','Peak') ;
          ShowCircularity(h.Ax_Gen3,h.NewShape3,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Circular_Skel_Clus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.BW,'Skeleton','Clustering') ;
          ShowCircularity(h.Ax_Ana,h.BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape1,'Skeleton','Clustering') ;
          ShowCircularity(h.Ax_Gen1,h.NewShape1,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape2,'Skeleton','Clustering') ;
          ShowCircularity(h.Ax_Gen2,h.NewShape2,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape3,'Skeleton','Clustering') ;
          ShowCircularity(h.Ax_Gen3,h.NewShape3,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Circular_Skel_Peak}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.BW,'Skeleton','Peak') ;
          ShowCircularity(h.Ax_Ana,h.BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape1,'Skeleton','Peak') ;
          ShowCircularity(h.Ax_Gen1,h.NewShape1,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape2,'Skeleton','Peak') ;
          ShowCircularity(h.Ax_Gen2,h.NewShape2,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xmic,ymic,rmic,xmcc,ymcc,rmcc,CircFactor] = Circularity(h.NewShape3,'Skeleton','Peak') ;
          ShowCircularity(h.Ax_Gen3,h.NewShape3,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,CircFactor,3) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.MZC}
      Processing('p',h.Process)
      if h.FlagSteps
          if h.Shape == 0
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.BW,'Centroid','Peak') ;
            ShowMZC(h.Ax_Ana,h.BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          elseif h.Shape == 1 && h.FlagGeneShape == 1
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.NewShape1,'Centroid','Peak') ;
            ShowMZC(h.Ax_Gen1,h.NewShape1,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          elseif h.Shape == 2 && h.FlagGeneShape == 1
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.NewShape2,'Centroid','Peak') ;
            ShowMZC(h.Ax_Gen2,h.NewShape2,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          elseif h.Shape == 3 && h.FlagGeneShape == 1
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.NewShape3,'Centroid','Peak') ;
            ShowMZC(h.Ax_Gen3,h.NewShape3,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          end
        end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.LSC}
      Processing('p',h.Process)
      if h.FlagSteps
          if h.Shape == 0
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.BW,'Centroid','Peak') ; ;
            ShowLSC(h.Ax_Ana,h.BW,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          elseif h.Shape == 1 && h.FlagGeneShape == 1
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.NewShape1,'Centroid','Peak') ; ;
            ShowLSC(h.Ax_Gen1,h.NewShape1,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          elseif h.Shape == 2 && h.FlagGeneShape == 1
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.NewShape1,'Centroid','Peak') ; ;
            ShowLSC(h.Ax_Gen2,h.NewShape2,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          elseif h.Shape == 3 && h.FlagGeneShape == 1
            [xmic , ymic , rmic , xmcc , ymcc , rmcc , CircFactor] = Circularity(h.NewShape1,'Centroid','Peak') ; ;
            ShowLSC(h.Ax_Gen3,h.NewShape3,xmic,ymic,rmic,xmcc,ymcc,rmcc,h.theta)
          end
        end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Compact}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,Perimeter,Area,CompactFactor] = Compactness(h.BW) ;
          ShowCompactness(h.Ax_Ana,h.BW,x,y)
          h.Table = Table_Update(h.Table,h.Shape,CompactFactor,4) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,Perimeter,Area,CompactFactor] = Compactness(h.NewShape1) ;
          ShowCompactness(h.Ax_Gen1,h.NewShape1,x,y)
          h.Table = Table_Update(h.Table,h.Shape,CompactFactor,4) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,Perimeter,Area,CompactFactor] = Compactness(h.NewShape2) ;
          ShowCompactness(h.Ax_Gen2,h.NewShape2,x,y)
          h.Table = Table_Update(h.Table,h.Shape,CompactFactor,4) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,Perimeter,Area,CompactFactor] = Compactness(h.NewShape3) ;
          ShowCompactness(h.Ax_Gen3,h.NewShape3,x,y)
          h.Table = Table_Update(h.Table,h.Shape,CompactFactor,4) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Convexity}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,ShapePerimeter,ConvexPerimeter,ConvexBoundary,ConvexFactor] = Convexity(h.BW) ;
          ShowConvexity(h.Ax_Ana,h.BW,x,y,ConvexBoundary)
          h.Table = Table_Update(h.Table,h.Shape,ConvexFactor,5) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,ShapePerimeter,ConvexPerimeter,ConvexBoundary,ConvexFactor] = Convexity(h.NewShape1) ;
          ShowConvexity(h.Ax_Gen1,h.NewShape1,x,y,ConvexBoundary)
          h.Table = Table_Update(h.Table,h.Shape,ConvexFactor,5) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,ShapePerimeter,ConvexPerimeter,ConvexBoundary,ConvexFactor] = Convexity(h.NewShape2) ;
          ShowConvexity(h.Ax_Gen2,h.NewShape2,x,y,ConvexBoundary)
          h.Table = Table_Update(h.Table,h.Shape,ConvexFactor,5) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,ShapePerimeter,ConvexPerimeter,ConvexBoundary,ConvexFactor] = Convexity(h.NewShape3) ;
          ShowConvexity(h.Ax_Gen3,h.NewShape3,x,y,ConvexBoundary)
          h.Table = Table_Update(h.Table,h.Shape,ConvexFactor,5) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Eccentricity}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [a,b,c,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory,EccentFactor] = Eccentricity(h.BW) ;
          ShowEccentricity(h.Ax_Ana,h.BW,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory)
          h.Table = Table_Update(h.Table,h.Shape,EccentFactor,6) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [a,b,c,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory,EccentFactor] = Eccentricity(h.NewShape1) ;
          ShowEccentricity(h.Ax_Gen1,h.NewShape1,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory)
          h.Table = Table_Update(h.Table,h.Shape,EccentFactor,6) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [a,b,c,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory,EccentFactor] = Eccentricity(h.NewShape2) ;
          ShowEccentricity(h.Ax_Gen2,h.NewShape2,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory)
          h.Table = Table_Update(h.Table,h.Shape,EccentFactor,6) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [a,b,c,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory,EccentFactor] = Eccentricity(h.NewShape3) ;
          ShowEccentricity(h.Ax_Gen3,h.NewShape3,Elli_Bound,X1,X2,Y1,Y2,x1,x2,y1,y2,xc,yc,onMajorx,onMajory,onMinorx,onMinory)
          h.Table = Table_Update(h.Table,h.Shape,EccentFactor,6) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Elliptic}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,ellix,elliy,elliperi,shapeperi,EllipFactor] = EllipticSmoothness(h.BW) ;
          ShowElliptic(h.Ax_Ana,h.BW,x,y,ellix,elliy)
          h.Table = Table_Update(h.Table,h.Shape,EllipFactor,7) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,ellix,elliy,elliperi,shapeperi,EllipFactor] = EllipticSmoothness(h.NewShape1) ;
          ShowElliptic(h.Ax_Gen1,h.NewShape1,x,y,ellix,elliy)
          h.Table = Table_Update(h.Table,h.Shape,EllipFactor,7) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,ellix,elliy,elliperi,shapeperi,EllipFactor] = EllipticSmoothness(h.NewShape2) ;
          ShowElliptic(h.Ax_Gen2,h.NewShape2,x,y,ellix,elliy)
          h.Table = Table_Update(h.Table,h.Shape,EllipFactor,7) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,ellix,elliy,elliperi,shapeperi,EllipFactor] = EllipticSmoothness(h.NewShape3) ;
          ShowElliptic(h.Ax_Gen3,h.NewShape3,x,y,ellix,elliy)
          h.Table = Table_Update(h.Table,h.Shape,EllipFactor,7) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Elongation}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x1,y1,x2,y2,x3,y3,x4,y4,major,minor,ElongFactor] = Elongation(h.BW) ;
          ShowElongation(h.Ax_Ana,h.BW,x1,x2,x3,x4,y1,y2,y3,y4)
          h.Table = Table_Update(h.Table,h.Shape,ElongFactor,8) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x1,y1,x2,y2,x3,y3,x4,y4,major,minor,ElongFactor] = Elongation(h.NewShape1) ;
          ShowElongation(h.Ax_Gen1,h.NewShape1,x1,x2,x3,x4,y1,y2,y3,y4)
          h.Table = Table_Update(h.Table,h.Shape,ElongFactor,8) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x1,y1,x2,y2,x3,y3,x4,y4,major,minor,ElongFactor] = Elongation(h.NewShape2) ;
          ShowElongation(h.Ax_Gen2,h.NewShape2,x1,x2,x3,x4,y1,y2,y3,y4)
          h.Table = Table_Update(h.Table,h.Shape,ElongFactor,8) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x1,y1,x2,y2,x3,y3,x4,y4,major,minor,ElongFactor] = Elongation(h.NewShape3) ;
          ShowElongation(h.Ax_Gen3,h.NewShape3,x1,x2,x3,x4,y1,y2,y3,y4)
          h.Table = Table_Update(h.Table,h.Shape,ElongFactor,8) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Regularity}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,xc,yc,sp,cp,RegularFactor] = Regularity(h.BW) ;
          ShowRegularity(h.Ax_Ana,h.BW,xc,yc,x,y)
          h.Table = Table_Update(h.Table,h.Shape,RegularFactor,9) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,xc,yc,sp,cp,RegularFactor] = Regularity(h.NewShape1) ;
          ShowRegularity(h.Ax_Gen1,h.NewShape1,xc,yc,x,y)
          h.Table = Table_Update(h.Table,h.Shape,RegularFactor,9) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,xc,yc,sp,cp,RegularFactor] = Regularity(h.NewShape2) ;
          ShowRegularity(h.Ax_Gen2,h.NewShape2,xc,yc,x,y)
          h.Table = Table_Update(h.Table,h.Shape,RegularFactor,9) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,xc,yc,sp,cp,RegularFactor] = Regularity(h.NewShape3) ;
          ShowRegularity(h.Ax_Gen3,h.NewShape3,xc,yc,x,y)
          h.Table = Table_Update(h.Table,h.Shape,RegularFactor,9) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Round_Cent_WClus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.BW,'WithClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Ana,h.BW,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape1,'WithClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Gen1,h.NewShape1,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape2,'WithClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Gen2,h.NewShape2,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape3,'WithClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Gen3,h.NewShape3,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Round_Cent_WOClus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.BW,'WithoutClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Ana,h.BW,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape1,'WithoutClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Gen1,h.NewShape1,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape2,'WithoutClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Gen2,h.NewShape2,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape3,'WithoutClust',0.98,'Centroid',false) ;
          ShowRoundness(h.Ax_Gen3,h.NewShape3,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Round_Skel_WClus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.BW,'WithClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Ana,h.BW,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape1,'WithClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Gen1,h.NewShape1,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape2,'WithClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Gen2,h.NewShape2,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape3,'WithClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Gen3,h.NewShape3,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Round_Skel_WOClus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.BW,'WithoutClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Ana,h.BW,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape1,'WithoutClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Gen1,h.NewShape1,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape2,'WithoutClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Gen2,h.NewShape2,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xcent,ycent,r,XMIC,YMIC,RMIC,RoundFactor]= Roundness(h.NewShape3,'WithoutClust',0.98,'Skeleton',false) ;
          ShowRoundness(h.Ax_Gen3,h.NewShape3,xcent,ycent,r,XMIC,YMIC,RMIC,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,RoundFactor,10) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.ShapeFactor_C_Outside}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.BW) ;
          ShowShapeFactor(h.Ax_Ana,h.BW,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape1) ;
          ShowShapeFactor(h.Ax_Gen1,h.NewShape1,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape2) ;
          ShowShapeFactor(h.Ax_Gen2,h.NewShape2,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape3) ;
          ShowShapeFactor(h.Ax_Gen3,h.NewShape3,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.ShapeFactor_C_Inside}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.BW,'Inside') ;
          ShowShapeFactor(h.Ax_Ana,h.BW,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape1,'Inside') ;
          ShowShapeFactor(h.Ax_Gen1,h.NewShape1,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape2,'Inside') ;
          ShowShapeFactor(h.Ax_Gen2,h.NewShape2,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape3,'Inside') ;
          ShowShapeFactor(h.Ax_Gen3,h.NewShape3,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_c,11) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.ShapeFactor_S_Outside}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.BW) ;
          ShowShapeFactor(h.Ax_Ana,h.BW,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape1) ;
          ShowShapeFactor(h.Ax_Gen1,h.NewShape1,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape2) ;
          ShowShapeFactor(h.Ax_Gen2,h.NewShape2,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape3) ;
          ShowShapeFactor(h.Ax_Gen3,h.NewShape3,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.ShapeFactor_S_Inside}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.BW,'Inside') ;
          ShowShapeFactor(h.Ax_Ana,h.BW,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape1,'Inside') ;
          ShowShapeFactor(h.Ax_Gen1,h.NewShape1,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape2,'Inside') ;
          ShowShapeFactor(h.Ax_Gen2,h.NewShape2,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,xc,yc,X1,X2,Y1,Y2,l,SF_c,SF_r] = ShapeFactor(h.NewShape3,'Inside') ;
          ShowShapeFactor(h.Ax_Gen3,h.NewShape3,x,y,xc,yc,X1,X2,Y1,Y2)
          h.Table = Table_Update(h.Table,h.Shape,SF_r,11) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Solidity}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [x,y,xc,yc,shapeA,convexA,SolidFactor] = Solidity(h.BW) ;
          ShowSolidity(h.Ax_Ana,h.BW,x,y,xc,yc)
          h.Table = Table_Update(h.Table,h.Shape,SolidFactor,12) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [x,y,xc,yc,shapeA,convexA,SolidFactor] = Solidity(h.NewShape1) ;
          ShowSolidity(h.Ax_Gen1,h.NewShape1,x,y,xc,yc)
          h.Table = Table_Update(h.Table,h.Shape,SolidFactor,12) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [x,y,xc,yc,shapeA,convexA,SolidFactor] = Solidity(h.NewShape2) ;
          ShowSolidity(h.Ax_Gen2,h.NewShape2,x,y,xc,yc)
          h.Table = Table_Update(h.Table,h.Shape,SolidFactor,12) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [x,y,xc,yc,shapeA,convexA,SolidFactor] = Solidity(h.NewShape3) ;
          ShowSolidity(h.Ax_Gen3,h.NewShape3,x,y,xc,yc)
          h.Table = Table_Update(h.Table,h.Shape,SolidFactor,12) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Sph_Without_Clus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.BW,'NoMIC','Clustering') ;
          ShowSphericity(h.Ax_Ana,h.BW,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape1,'NoMIC','Clustering') ;
          ShowSphericity(h.Ax_Gen1,h.NewShape1,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape2,'NoMIC','Clustering') ;
          ShowSphericity(h.Ax_Gen2,h.NewShape2,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape3,'NoMIC','Clustering') ;
          ShowSphericity(h.Ax_Gen3,h.NewShape3,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Sph_Without_Peak}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.BW,'NoMIC','Peak') ;
          ShowSphericity(h.Ax_Ana,h.BW,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape1,'NoMIC','Peak') ;
          ShowSphericity(h.Ax_Gen1,h.NewShape1,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape2,'NoMIC','Peak') ;
          ShowSphericity(h.Ax_Gen2,h.NewShape2,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape3,'NoMIC','Peak') ;
          ShowSphericity(h.Ax_Gen3,h.NewShape3,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Sph_WithCent_Clus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.BW,'Centroid','Clustering') ;
          ShowSphericity(h.Ax_Ana,h.BW,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape1,'Centroid','Clustering') ;
          ShowSphericity(h.Ax_Gen1,h.NewShape1,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape2,'Centroid','Clustering') ;
          ShowSphericity(h.Ax_Gen2,h.NewShape2,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape3,'Centroid','Clustering') ;
          ShowSphericity(h.Ax_Gen3,h.NewShape3,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Sph_WithCent_Peak}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.BW,'Centroid','Peak') ;
          ShowSphericity(h.Ax_Ana,h.BW,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape1,'Centroid','Peak') ;
          ShowSphericity(h.Ax_Gen1,h.NewShape1,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape2,'Centroid','Peak') ;
          ShowSphericity(h.Ax_Gen2,h.NewShape2,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape3,'Centroid','Peak') ;
          ShowSphericity(h.Ax_Gen3,h.NewShape3,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Sph_WithSkel_Clus}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.BW,'Skeleton','Clustering') ;
          ShowSphericity(h.Ax_Ana,h.BW,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape1,'Skeleton','Clustering') ;
          ShowSphericity(h.Ax_Gen1,h.NewShape1,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape2,'Skeleton','Clustering') ;
          ShowSphericity(h.Ax_Gen2,h.NewShape2,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape3,'Skeleton','Clustering') ;
          ShowSphericity(h.Ax_Gen3,h.NewShape3,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Sph_WithSkel_Peak}
      Processing('p',h.Process)
      if h.FlagSteps
        if h.Shape == 0
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.BW,'Skeleton','Peak') ;
          ShowSphericity(h.Ax_Ana,h.BW,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 1 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape1,'Skeleton','Peak') ;
          ShowSphericity(h.Ax_Gen1,h.NewShape1,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 2 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape2,'Skeleton','Peak') ;
          ShowSphericity(h.Ax_Gen2,h.NewShape2,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        elseif h.Shape == 3 && h.FlagGeneShape == 1
          [xc,yc,rmcc,A,r_equivalent,SpherFactor]= Sphericity(h.NewShape3,'Skeleton','Peak') ;
          ShowSphericity(h.Ax_Gen3,h.NewShape3,xc,yc,rmcc,h.theta)
          h.Table = Table_Update(h.Table,h.Shape,SpherFactor,13) ;
        end
      end

      guidata(gcf,h) % update handles
      Processing('r',h.Process)

    case {h.Run}
      Processing('p',h.Process)
      if h.FlagSteps
        try
          ModuleSelector ;
          uiwait ;
          h.GeneMethod = getappdata(0,'GeneMethod') ;
          ####### Processing
          % Get Boundaries
          [Bound,Lab,~] = bwboundaries(h.BW) ;
          X = Bound{1,1}(:,2) ;
          Y = Bound{1,1}(:,1) ;
          % Get Centroid
          S = regionprops(h.BW,'Centroid') ;
          % Get Distance
          for i = 1:length(X)
            D(i) = pdist([S.Centroid(1) S.Centroid(2) ; X(i) Y(i)]) ;
          end
          D = D' ;
          % Get Angle
          t = AngleFinder(S.Centroid(1),S.Centroid(2),X,Y) ;
          % Rearrange Data
          [X,Y,D,t] = ReArrange(X,Y,D,t) ;
          % Smooth Data
          D = Smoother(D,0.01) ; % 0.01
          % Normalize & Find threshold normalized data
          NewD = normalize(D,'2') ;
          W = mean(D./NewD) ;
          % Find Extremas
          [BaseLoc , MaxLoc , MaxExt , MinLoc , MinExt] = findextremas(NewD) ;
          % Extremas Indices
          Base = sort(unique([BaseLoc ; MinLoc])) ;
          Base = [1 ; Base ; length(NewD)] ;
          [index , K] = ExtremasToIndex(Base,MaxLoc,MaxExt,MinExt,NewD) ;
          % Lower and Upper Conditoions For Optimization
          xD = [1:length(NewD)]' ;
          [Lower , Upper , Guess] = BoundFinder(NewD,K,index) ;
          % CurveFit(lsqcurvefit) run time 01min:00sec
          options = optimset('Algorithm','lsqcurvefit','Display','off','MaxIter',50000,'TolFun',1e-14,'FinDiffRelStep',1e-14) ; % 'MaxFunctionEvaluations',5000,'FiniteDifferenceStepSize',1e-16,'OptimalityTolerance',1e-16,'StepTolerance',1e-16 ,'OptimalityTolerance',1e-16,'MaxFunctionEvaluations',100000
          tic
          Coeffs = lsqcurvefit(@ObjectiveFunc,Guess,xD,NewD,Lower,Upper,options) ;
          time = toc ; % sec
          % Reshape Coefficients Matrix
          AINumber = reshape(Coeffs,[K,3]) ;
          % Distance Calculation
          FitD = FObjectiveFunc(xD,AINumber) ;
          % Initialize Original Shape For Export
          Rotate = 0 ; % Degree
          [ScaledXShow,ScaledYShow,~,~] = BoundGenerator(S.Centroid(1),S.Centroid(2),FitD,W,t,Rotate) ;
          [~,~,NewBWShow] = XY2BWMatrix(ScaledXShow,ScaledYShow,size(h.BW)) ;
          % Initialize Generated Shapes For Export
          % Generate Random Index For New Shapes
          NewAINumber_index = FindIndex2NewAINumbers(K) ;
          % Select a random model to generate shapes
          [AINumber_S1 , AINumber_S2 , AINumber_S3] = ModelSelector(AINumber,NewAINumber_index,h.GeneMethod) ;
          Shape1 = FObjectiveFunc(xD,AINumber_S1) ;
          Shape2 = FObjectiveFunc(xD,AINumber_S2) ;
          Shape3 = FObjectiveFunc(xD,AINumber_S3) ;
          [ScaledX,ScaledY,~,~] = BoundGenerator(0,0,FitD,W,t,0) ;
          [NewScaledX,NewScaledY,NewBW] = XY2BWMatrix(ScaledX,ScaledY,size(h.BW)) ;
          [ScaledX_S1,ScaledY_S1,~,~] = BoundGenerator(0,0,Shape1,W,t,0) ;
          [NewScaledX_S1,NewScaledY_S1,NewShape1] = XY2BWMatrix(ScaledX_S1,ScaledY_S1,size(h.BW)) ;
          imshow(~NewShape1,'parent',h.Ax_Gen1)
          [ScaledX_S2,ScaledY_S2,~,~] = BoundGenerator(0,0,Shape2,W,t,0) ;
          [NewScaledX_S2,NewScaledY_S2,NewShape2] = XY2BWMatrix(ScaledX_S2,ScaledY_S2,size(h.BW)) ;
          imshow(~NewShape2,'parent',h.Ax_Gen2)
          [ScaledX_S3,ScaledY_S3,~,~] = BoundGenerator(0,0,Shape3,W,t,0) ;
          [NewScaledX_S3,NewScaledY_S3,NewShape3] = XY2BWMatrix(ScaledX_S3,ScaledY_S3,size(h.BW)) ;
          imshow(~NewShape3,'parent',h.Ax_Gen3)
          % Save Binary Image
          h.NewBW = NewBW ;
          h.NewBWShow = NewBWShow ;
          h.NewShape1 = NewShape1 ;
          h.NewShape2 = NewShape2 ;
          h.NewShape3 = NewShape3 ;
          % Save Structure
          Org.XC = 0 ; Org.YC = 0 ; Org.W = W ; Org.X = ScaledX ; Org.Y = ScaledY ; Org.D = FitD ; Org.Theta = t ;
          Sh1.XC = 0 ; Sh1.YC = 0 ; Sh1.W = W ; Sh1.X = ScaledX_S1 ; Sh1.Y = ScaledY_S1 ; Sh1.D = Shape1 ; Sh1.Theta = t ;
          Sh2.XC = 0 ; Sh2.YC = 0 ; Sh2.W = W ; Sh2.X = ScaledX_S2 ; Sh2.Y = ScaledY_S2 ; Sh2.D = Shape2 ; Sh2.Theta = t ;
          Sh3.XC = 0 ; Sh3.YC = 0 ; Sh3.W = W ; Sh3.X = ScaledX_S3 ; Sh3.Y = ScaledY_S3 ; Sh3.D = Shape3 ; Sh3.Theta = t ;
          h.Org = Org ; h.Sh1 = Sh1 ; h.Sh2 = Sh2 ; h.Sh3 = Sh3 ;
          h.FlagGeneShape = 1 ;
          Processing('r',h.Process)
          ####### End of Processing
        catch
          Processing('e',h.Process)
        end

      end

      guidata(gcf,h) % update handles

    case {h.SaveGeneratedImage}
      Processing('p',h.Process)
      if h.FlagGeneShape
        [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group"}) ;
        if (FileIndex) ~= 0
          F1 = getframe(h.Ax_Gen1);
          Image1 = frame2im(F1);
          F2 = getframe(h.Ax_Gen2);
          Image2 = frame2im(F2);
          F3 = getframe(h.Ax_Gen3);
          Image3 = frame2im(F3);
          imwrite(Image1, [FilePath , strsplit(FileName,'.'){1:end-1} , '1' , '.' , strsplit(FileName,'.'){end}])
          imwrite(Image2, [FilePath , strsplit(FileName,'.'){1:end-1} , '2' , '.' , strsplit(FileName,'.'){end}])
          imwrite(Image3, [FilePath , strsplit(FileName,'.'){1:end-1} , '3' , '.' , strsplit(FileName,'.'){end}])
        end
      end

      Processing('r',h.Process)

    case {h.SaveGeneratedJSON}
      Processing('p',h.Process)
      if h.FlagGeneShape
        [FileName, FilePath, FileIndex] = uiputfile({"*.json","JavaScript Object Notation"}) ;
        if (FileIndex) ~= 0
          SturctToJSON({h.Sh1;h.Sh2;h.Sh3},[FilePath FileName])
        end
      end

      Processing('r',h.Process)

    case {h.RunPM}

      if h.FlagGeneShape
        try
          Processing('p',h.Process)
          GetDataFromUser ;
          uiwait ;
          h.Dimension = getappdata(0,'Dimension') ;
          h.Sorting = getappdata(0,'Sorting') ;
          h.Packing = getappdata(0,'Packing') ;
          [R1,R5,R16,R50,R84,R95,R100,a,x,Fun,Category,xx,yy,PhiValue] = SortRequest(h.Sorting) ;
          %% Constructing a k-item distribution of grains
          k = 10 ;
          r = GenerateKthRadius(xx,yy,k) ;
          Dim = h.Dimension ;
          R_min = min(r) ;
          R_max = max(r) ;
          Bound = true ;
          if R_max.^2 < Dim(1)*Dim(2)
              [C,R] = RandomShapeGenerator(Dim,R_min,R_max,Bound) ;
          else
              Processing('e',h.Process)
          end
          %% Packing
          R = R-0.5*R*(h.Packing) ;
          Rotate = randi(360,size(R)) ; % Degree [0 359]
          DC = [(h.Org.D).' ; (h.Sh1.D).' ; (h.Sh2.D).' ; (h.Sh2.D).'] ;
          TC = [(h.Org.Theta).' ; (h.Sh1.Theta).' ; (h.Sh2.Theta).' ; (h.Sh2.Theta).'] ;
          INDI = randi(size(DC,1),[size(R,1) 1]) ;
          for i = 1:length(R)
            DD(i,:) = DC(INDI(i),:) ;
            T(i,:) = TC(INDI(i),:) ;
          end
          Coef = R./max(DD,[],2) ;
          R = DD.*Coef ;
          [X,Y] = ShapePainter(C(:,1),C(:,2),R,T,Rotate) ;
          Shapes.XC = C(:,1) ; Shapes.YC = C(:,2) ; Shapes.W = 1 ; Shapes.X = X ; Shapes.Y = Y ; Shapes.D = R ; Shapes.Theta = T+180+Rotate ;
          h.Shapes = Shapes ;
          h.FlagGenePM = 1 ;
          Processing('r',h.Process)
        catch
          Processing('e',h.Process)
        end
      end

      guidata(gcf,h) % update handles

    case {h.Save2DPMImage}
      Processing('p',h.Process)
      if h.FlagGenePM
        [FileName, FilePath, FileIndex] = uiputfile({"*.jpeg;*.jpg;*.png", "Supported Picture Formats";"*.png", "Portable Network Graphics";"*.jpeg ; *.jpg", "Joint Photographic Experts Group"}) ;
        if (FileIndex) ~= 0
          figure('name',"Saved Image",'NumberTitle','off','resize','off',"toolbar", "none",'uicontextmenu',[],'menubar','none') ; Ax = axes("Units",'Normalized',"Position",[0 0 1 1]) ;
          for i = 1:size(h.Shapes.X ,1)
            hold on ; patch(Ax,h.Shapes.X(i,:),h.Shapes.Y(i,:),"k",'EdgeColor',"none")
          end
          set(Ax,'XColor','none') ; set(Ax,'YColor','none')
          axis([0 h.Dimension(1) 0 h.Dimension(2)])
          F = getframe(Ax);
          Image = frame2im(F);
          imwrite(Image, [FilePath ,  FileName])
        end
      end

      Processing('r',h.Process)

    case {h.Save2DPMJSON}
      Processing('p',h.Process)
      if h.FlagGenePM
        [FileName, FilePath, FileIndex] = uiputfile({"*.json","JavaScript Object Notation"}) ;
        if (FileIndex) ~= 0
          Cell = {} ;
          for i = 1:size(h.Shapes.X ,1)
            Struc.XC = h.Shapes.XC(i,1) ; Struc.YC = h.Shapes.YC(i,1) ; Struc.W = 1 ; Struc.X = h.Shapes.X(i,:) ; Struc.Y = h.Shapes.Y(i,:) ; Struc.D = h.Shapes.D(i,:) ; Struc.Theta = h.Shapes.Theta(i,:) ;
            Cell{i,1} = Struc ;
          end
          SturctToJSON(Cell,[FilePath FileName])
        end
      end

      Processing('r',h.Process)

  end % Switch
end

##
root.Fig = figure("toolbar", "none",'uicontextmenu',[],'menubar','none','name',"Random Porous Medium Generator 2D",'NumberTitle','off','resize','off','units','normalized',"Position", [2.1962e-03 6.2500e-02 9.9414e-01 8.6979e-01],"CloseRequestFcn",'exit') ; % for .exe
root.F = uimenu("label", "&File", "accelerator", "f",'callback',@Update_UI);
root.E = uimenu("label", "&Edit", "accelerator", "e");
root.A = uimenu("label", "&Analysis", "accelerator", "a");
root.G = uimenu("label", "&Generator", "accelerator", "g");
root.H = uimenu("label", "&Help", "accelerator", "h");

# Subs
root.Open = uimenu(root.F, "label", "&Open Image", "accelerator", "O", "callback", @Update_UI);
uimenu(root.F, "label", "E&xit", "accelerator", "X", "callback", 'close(root.Fig) ; exit'); % for .exe

root.SaveImgAn = uimenu(root.E, "label", "&Save Last Analysis as Image", "accelerator", "S","callback", @Update_UI);
root.SaveTable = uimenu(root.E, "label", "Save &Table as Sheet(.csv)", "accelerator", "T","callback", @Update_UI);
root.SaveImgRe = uimenu(root.E, "label", "Save Reconstructed Shape as Image", "accelerator", "","callback", @Update_UI);
root.SaveImgJSON = uimenu(root.E, "label", "Save Shape as JSON", "accelerator", "","callback", @Update_UI);

root.AngulMenu = uimenu(root.A, "label", "Angularity", "accelerator", "","callback", "");
root.Angularity_Cent = uimenu(root.AngulMenu, "label", "With Centroid Method", "accelerator", "","callback", @Update_UI);
root.Angularity_Skel = uimenu(root.AngulMenu, "label", "With Skeleton Method", "accelerator", "","callback", @Update_UI);

root.AspectMenu = uimenu(root.A, "label", "Aspect Ratio", "accelerator", "","callback", "");
root.AspectRatio_Outside = uimenu(root.AspectMenu, "label", "Axes anywhere", "accelerator", "","callback", @Update_UI);
root.AspectRatio_Inside  = uimenu(root.AspectMenu, "label", "Axes inside", "accelerator", "","callback", @Update_UI);

root.Circular = uimenu(root.A, "label", "Circularity", "accelerator", "","callback", "");
root.Circular_Cent = uimenu(root.Circular, "label", "MIC Centre With Centroid", "accelerator", "","callback", "");
root.Circular_Cent_Clus = uimenu(root.Circular_Cent, "label", "MCC With Distance Clustering", "accelerator", "","callback", @Update_UI);
root.Circular_Cent_Peak = uimenu(root.Circular_Cent, "label", "MCC Centre With Distance Peak", "accelerator", "","callback", @Update_UI);

root.Circular_Skel = uimenu(root.Circular, "label", "MIC Centre With Skeleton", "accelerator", "","callback", "");
root.Circular_Skel_Clus = uimenu(root.Circular_Skel, "label", "MCC With Distance Clustering", "accelerator", "","callback", @Update_UI);
root.Circular_Skel_Peak = uimenu(root.Circular_Skel, "label", "MCC Centre With Distance Peak", "accelerator", "","callback", @Update_UI);

root.MZC = uimenu(root.Circular, "label", "Minimum Zone Circle(MZC)", "accelerator", "","callback", @Update_UI);
root.LSC = uimenu(root.Circular, "label", "Least Squares Circle(LSC)", "accelerator", "","callback", @Update_UI);

root.Compact = uimenu(root.A, "label", "Compactness", "accelerator", "","callback", @Update_UI);

root.Convexity = uimenu(root.A, "label", "Convexity", "accelerator", "","callback", @Update_UI);

root.Eccentricity = uimenu(root.A, "label", "Eccentricity", "accelerator", "","callback", @Update_UI);

root.Elliptic = uimenu(root.A, "label", "Elliptic Smoothness", "accelerator", "","callback", @Update_UI);

root.Elongation = uimenu(root.A, "label", "Elongation", "accelerator", "","callback", @Update_UI);

root.Regularity = uimenu(root.A, "label", "Regularity", "accelerator", "","callback", @Update_UI);

root.Roundness = uimenu(root.A, "label", "Roundness", "accelerator", "","callback", "");
root.Round_Cent = uimenu(root.Roundness, "label", "MIC Centre With Centroid", "accelerator", "","callback", "");
root.Round_Cent_WClus = uimenu(root.Round_Cent, "label", "Corners with Cluster", "accelerator", "","callback", @Update_UI);
root.Round_Cent_WOClus = uimenu(root.Round_Cent, "label", "Corners without Cluster", "accelerator", "","callback", @Update_UI);

root.Round_Skel = uimenu(root.Roundness, "label", "MIC Centre With Skeleton", "accelerator", "","callback", "");
root.Round_Skel_WClus = uimenu(root.Round_Skel, "label", "Corners with Cluster", "accelerator", "","callback", @Update_UI);
root.Round_Skel_WOClus = uimenu(root.Round_Skel, "label", "Corners without Cluster", "accelerator", "","callback", @Update_UI);

root.ShapeMenu = uimenu(root.A, "label", "ShapeFactors", "accelerator", "","callback", "");
root.ShapeFactor_Circ = uimenu(root.ShapeMenu, "label", "Circularity", "accelerator", "","callback","");
root.ShapeFactor_C_Outside = uimenu(root.ShapeFactor_Circ, "label", "Major Axis anywhere", "accelerator", "","callback", @Update_UI);
root.ShapeFactor_C_Inside = uimenu(root.ShapeFactor_Circ, "label", "Major Axis inside", "accelerator", "","callback", @Update_UI);

root.ShapeFactor_Sphe = uimenu(root.ShapeMenu, "label", "Sphericity", "accelerator", "","callback","");
root.ShapeFactor_S_Outside = uimenu(root.ShapeFactor_Sphe, "label", "Major Axis anywhere", "accelerator", "","callback", @Update_UI);
root.ShapeFactor_S_Inside = uimenu(root.ShapeFactor_Sphe, "label", "Major Axis inside", "accelerator", "","callback", @Update_UI);

root.Solidity = uimenu(root.A, "label", "Solidity", "accelerator", "","callback", @Update_UI);

root.Sphericity = uimenu(root.A, "label", "Sphericity", "accelerator", "","callback", "");
root.Without = uimenu(root.Sphericity, "label", "Without MIC Centre", "accelerator", "","callback", "");
root.Sph_Without_Clus = uimenu(root.Without, "label", "MCC with Distance Clustering", "accelerator", "","callback", @Update_UI);
root.Sph_Without_Peak = uimenu(root.Without, "label", "MCC with Distance Peaks", "accelerator", "","callback", @Update_UI);

root.WithCent = uimenu(root.Sphericity, "label", "MIC Centre With Centroid", "accelerator", "","callback", "");
root.Sph_WithCent_Clus = uimenu(root.WithCent, "label", "MCC with Distance Clustering", "accelerator", "","callback", @Update_UI);
root.Sph_WithCent_Peak = uimenu(root.WithCent, "label", "MCC with Distance Peaks", "accelerator", "","callback", @Update_UI);

root.WithSkel = uimenu(root.Sphericity, "label", "MIC Centre With Skeleton", "accelerator", "","callback", "");
root.Sph_WithSkel_Clus = uimenu(root.WithSkel, "label", "MCC with Distance Clustering", "accelerator", "","callback", @Update_UI);
root.Sph_WithSkel_Peak = uimenu(root.WithSkel, "label", "MCC with Distance Peaks", "accelerator", "","callback", @Update_UI);

root.Run = uimenu(root.G, "label", "&Run to Generate New Shapes Like Original Shape", "accelerator", "R","callback", @Update_UI);
root.SaveGeneratedImage = uimenu(root.G, "label", "Save Generated Shapes as Image", "accelerator", "","callback", @Update_UI);
root.SaveGeneratedJSON = uimenu(root.G, "label", "Save Generated Shapes as JSON", "accelerator", "","callback", @Update_UI); % Octave 7.3
root.RunPM = uimenu(root.G, "label", "Run to Generate &Porous Medium With Generated Shapes", "accelerator", "P","callback", @Update_UI);
root.Save2DPMImage = uimenu(root.G, "label", "Save Porous Medium as Image", "accelerator", "","callback", @Update_UI);
root.Save2DPMJSON = uimenu(root.G, "label", "Save Porous Medium as JSON", "accelerator", "","callback", @Update_UI); % Octave 7.3


uimenu(root.H, "label", "&Documentation", "accelerator", "D","callback", "system(['start ./Reference/Manual.pdf']) ;");
uimenu(root.H, "label", "About Me", "accelerator", "","callback", "web('https://www.linkedin.com/in/seyed-mousa-sajadi-8284b1124/','-new')");

root.Ax_Org = axes("position", [0.015 0.510 0.315 0.425],'box','on','xtick',[],'ytick',[],'colormap',colormap('gray'));
root.Ax_Ana = axes("position", [0.015 0.020 0.315 0.425],'box','on','xtick',[],'ytick',[],'colormap',colormap('gray'));
root.txt_Org = uicontrol (root.Fig, "style", "text", "string","Original Image","units","normalized","position", [0.112 0.94 0.1 0.03],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"));
root.txt_Ana = uicontrol (root.Fig, "style", "text", "string","Latest Image Analysis","units","normalized","position", [0.112 0.455 0.1 0.03],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"));

root.P = uipanel ("title", "Analysis of Parameters for Original/Generated Shape", "position", [0.340 0.85 0.3195 0.10],'backgroundcolor',get(root.Fig,"Color"));
root.Radio0 = uicontrol (root.P, "style", "radiobutton", "string","Original   ","units","normalized","position", [0.02 0.35 0.317 0.35],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"),'value' , 1,'callback',@Update_UI);
root.Radio1 = uicontrol (root.P, "style", "radiobutton", "string","Generated 1","units","normalized","position", [0.25 0.35 0.317 0.35],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);
root.Radio2 = uicontrol (root.P, "style", "radiobutton", "string","Generated 2","units","normalized","position", [0.50 0.35 0.317 0.35],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);
root.Radio3 = uicontrol (root.P, "style", "radiobutton", "string","Generated 3","units","normalized","position", [0.75 0.35 0.317 0.35],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"),'value' , 0,'callback',@Update_UI);

root.TableData = nan(13,4);
root.Row_Names = { "Angularity", "Aspect Ratio", "Circularity" , "Compactness" , "Convexity" , "Eccentricity" , "Elliptic Smoothness" , "Elongation" , "Regularity" , "Roundness" , "ShapeFactor" , "Solidity" , "Sphericity"};
root.Col_Names = { "Original", "Generated 1", "Generated 2" , "Generated 3"};
root.Table = uitable (root.Fig, "Data", root.TableData, "RowName", root.Row_Names, "ColumnName", root.Col_Names,"units","normalized","Position",[0.340 0.365 0.3195 0.453],"fontsize" ,8);

root.Ax_Gen1 = axes("position", [0.67 0.660 0.315 0.28],'box','on','xtick',[],'ytick',[],'colormap',colormap('gray'));
root.Ax_Gen2 = axes("position", [0.67 0.340 0.315 0.28],'box','on','xtick',[],'ytick',[],'colormap',colormap('gray'));
root.Ax_Gen3 = axes("position", [0.67 0.020 0.315 0.28],'box','on','xtick',[],'ytick',[],'colormap',colormap('gray'));
root.txt_Gen1 = uicontrol (root.Fig, "style", "text", "string","Generated Shape 1","units","normalized","position", [0.79 0.945 0.1 0.03],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"));
root.txt_Gen2 = uicontrol (root.Fig, "style", "text", "string","Generated Shape 2","units","normalized","position", [0.79 0.625 0.1 0.03],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"));
root.txt_Gen3 = uicontrol (root.Fig, "style", "text", "string","Generated Shape 3","units","normalized","position", [0.79 0.305 0.1 0.03],'fontsize',8,'backgroundcolor',get(root.Fig,"Color"));
root.Process = uicontrol (root.Fig, "style", "text", "string"," Ready","units","normalized","position", [0.00 0.00 0.055 0.018],'fontsize',6,'backgroundcolor',[0.785 1 0.785],"horizontalalignment",'left');


guidata(gcf,root) ;
Update_UI(gcf,true)
pause


# cd D:\Full_Codes\Octave\RPMG2D
# cls ; octave .\RPMG2D.m
