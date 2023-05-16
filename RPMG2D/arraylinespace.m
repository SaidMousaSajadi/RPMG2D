function A = arraylinespace(A,Space,maxsize)
  a = [] ;
  for i = 1:size(A,1)
    Back = A(i)-Space ;
    Front = A(i)+Space ;
    Array = Back:1:Front ;
    Array(Array<1) = Array(Array<1)+maxsize ;
    Array(Array>maxsize) = Array(Array>maxsize)-maxsize ;
    a = [a Array] ;
  end
  A = sort(unique(a))' ;
end
