function K = evalcluster(out,S)
  % The average distance of the members to the center of the cluster (S)
  XCON = out(:,1) ;
  YCON = out(:,2) ;
  KList = [2:length(XCON)] ;
  for k = 1:length(KList)-1
    index = kmeans(out,KList(k)) ;
    for j = 1:KList(k)
      Variable(index==j,1:2) = repmat([mean(XCON(index==j)) , mean(YCON(index==j))],size(out(index==j,:),1),1) ;
    end
    Variable(:,3) = sqrt(sum((out(:,1:2) - Variable(:,1:2)).^2,2)) ;
    KListMean(k) = mean(Variable(:,3)) ;
    if KListMean(k) < S
      K = KList(k) ;
      break
    end
  end

end
