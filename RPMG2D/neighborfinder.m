function C = neighborfinder(X,Y,XCONi,YCONi,NeighborPixel) ;
  OO = [X,Y] ;
  AA = [XCONi YCONi] ;
  A = [find(all(OO-AA==[0 0],2))] ;
  C = arraylinespace(A,NeighborPixel,size(X,1)) ;
end
