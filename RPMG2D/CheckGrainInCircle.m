function varargout = CheckGrainInCircle(x,y,xc,yc,r,ThresholdFactor)

  Factor = nnz(r > sqrt((x-xc).^2 + (y-yc).^2))/length(x) ;

  varargout{1} = true ;
  varargout{2} = Factor ;
  if Factor < ThresholdFactor
      varargout{1} = false ;
  end

end % function
