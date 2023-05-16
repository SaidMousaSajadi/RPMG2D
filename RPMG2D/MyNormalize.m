function Xnew = MyNormalize(X,Method)
  switch lower(Method)
    case {'range'}
      if length(unique(X)) > 1
        a = 0 ; b = 1 ;
        Xnew = a + ((X-min(X(:)))/(max(X(:))-min(X(:))))*(b-a) ;
      else
        Xnew = ones(size(X)) ;
      end
##  case {'sort'}
##    [~,ii]=sort(X,'descend') ;
##    for i = 1:length(X)
##      S(ii(i)) = length(X)+1-i ;
##    endfor
##    a = 0 ; b = 1 ;
##    Xnew = a + ((S-min(S(:)))/(max(S(:))-min(S(:))))*(b-a) ;
  endswitch
end
