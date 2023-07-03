function mycsvwrite(FullPath,Row,Col,Data)
  Col = fliplr(Col) ; Col(1,size(Col,2)+1) = "   " ; Col = fliplr(Col) ;
  GetVer = [ispc ismac isunix] ;
  if GetVer == [1 0 0] % MSWindows
    NewLine = "\r" ;
  elseif GetVer == [0 1 0] % MacOS
    NewLine = "\r" ;
  elseif GetVer == [0 0 1] % Linux
    NewLine = "\r" ;
  else
    NewLine = "\n" ;
    disp("I`dont know")
  end
  CSVID = fopen(FullPath,'w+') ;
  % Header
  for i = 1:size(Col,2)
    csvwrite(CSVID, Col{i}, 'delimiter', ',',  'precision', '%s',"append",'on') ;
  end
  csvwrite(CSVID, "", 'delimiter', ',',  'precision', NewLine,"append",'on') ;
  % Rows
  for i = 1:size(Row,2)
    csvwrite(CSVID, Row{i}, 'delimiter', ',',  'precision', '%s',"append",'on') ;
    for j = 1:size(Data,2)
      csvwrite(CSVID, num2str(Data(i,j),"%3.5f"), 'delimiter', ',',  'precision', '%s',"append",'on') ;
    end
    csvwrite(CSVID, "", 'delimiter', ',',  'precision', NewLine,"append",'on') ;
  end
  fclose(CSVID) ;
##  fclose all ;
  end
