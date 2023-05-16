function L = CheckCenterInGrain(bw,xc,yc)
  try 
    Cent_Pixel = bw(round(yc),round(xc)) ;
    if Cent_Pixel == 1
      L = true ;
    else
      L = false ;
    end
  catch
    L = false ;
  end_try_catch
  
end % function