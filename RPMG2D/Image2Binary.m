function [BW,M] = Image2Binary(IM,M)
  pkg load image ;
  if ndims(IM) == 2 % 2D
    if islogical(IM) % Binary
      BW = IM(:,:,1) ;
    else % No-Binary
      if isempty(M) % No-ColorMap
        G = double(IM)/255 ;
        BW = im2bw(G,"moments") ;
      else % ColorMap
        IM = ind2rgb(IM,M) ;
        G = rgb2gray(IM) ;
        BW = im2bw(G,"moments") ;
      end
    end
  elseif ndims(IM) == 3 & size(IM,3) == 3 % 2D & 3 Channel
    if islogical(IM) % Binary
      BW = IM(:,:,1) ;
    else % No-Binary
      if isempty(M) % No-ColorMap
        G = rgb2gray(IM) ;
        BW = im2bw(G,"moments") ;
      else % ColorMap
        IM = ind2rgb(IM,M) ;
        G = rgb2gray(IM) ;
        BW = im2bw(G,"moments") ;
      end
    end
  end
end
