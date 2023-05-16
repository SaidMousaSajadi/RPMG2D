function varargout = CheckCircleInGrain(x,y,xc,yc,r,ThresholdFactor)
  % x , y is boundary of grain
  Factor = nnz(r < sqrt((x-xc).^2 + (y-yc).^2))/length(x) ;
  varargout{1} = false ;
  varargout{2} = Factor ;
  if Factor >= ThresholdFactor
      varargout{1} = true ;
  end
end % function
