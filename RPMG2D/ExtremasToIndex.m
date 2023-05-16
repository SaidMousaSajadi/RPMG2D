function [index , K] = ExtremasToIndex(Base,MaxLoc,MaxExt,MinExt,NewD)
  index = zeros(size(NewD)) ;
  K = 0 ;
  for i = 1:length(Base)-1
    if intersect(Base(i):Base(i+1),MaxLoc)
      K = K+1 ;
      if i == 1
        Flag = 1 ;
      elseif i == length(Base)-1
        Flag = 3 ;
      else
        Flag = 2 ;
      end
    elseif intersect(Base(i):Base(i+1),MaxExt)
      K = K+1 ;
      Flag = 4 ;
    elseif intersect(Base(i):Base(i+1),MaxLoc) & intersect(Base(i):Base(i+1),MaxExt)
      K = K+1 ;
      if i == 1
        Flag = 1 ;
      elseif i == length(Base)-1
        Flag = 3 ;
      else
        Flag = 2 ;
      end
    elseif intersect(Base(i):Base(i+1),MinExt)
      K = K+1 ;
      Flag = 5 ;
    elseif intersect(Base(i):Base(i+1),MaxLoc) & intersect(Base(i):Base(i+1),MinExt)
      K = K+1 ;
      if i == 1
        Flag = 1 ;
      elseif i == length(Base)-1
        Flag = 3 ;
      else
        Flag = 2 ;
      end
    else
      Flag = 0 ;
    end
    switch Flag
      case 1
        index(Base(i):Base(i+1)) = K ;
      case 2
        index(Base(i):Base(i+1)) = K ;
      case 3
        index(Base(end-1):Base(end)) = K ;
      case 4
        index(Base(i):intersect(Base(i):Base(i+1),MaxExt)) = K ;
      case 5
        index(intersect(Base(i):Base(i+1),MinExt):Base(i+1)) = K ;
      otherwise
      end
    end
endfunction

