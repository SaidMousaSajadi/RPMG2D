function C = neighborsfinder(X,Y,XCON,YCON,NeighborPixel)
  OO = [X,Y] ;
  AA = [XCON YCON] ;
  A = [] ;
  for i = 1:size(AA,1)
    A = [A ; find(all(OO-AA(i,:)==[0 0],2))] ;
  end
  C = arraylinespace(A,NeighborPixel,size(X,1)) ;
end
