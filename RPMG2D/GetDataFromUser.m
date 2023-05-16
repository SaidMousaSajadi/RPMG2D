function GetDataFromUser() ;
  function GetData(obj,init = false)
    h = guidata(obj) ;
    switch (gcbo)
      case {h.Gene}
        setappdata(0,'Dimension',str2num(get(h.Dime,'string'))) ;
        setappdata(0,'Sorting',get(h.Sort,'Value')) ;
        setappdata(0,'Packing',str2num(get(h.Pack,'string'))) ;
        close(h.F)
    endswitch
  end

  R.F = figure("toolbar", "none",'uicontextmenu',[],'menubar','none','name',"Dimensions, sorting, packing 2D model",'NumberTitle','off','resize','off','units','normalized',"Position", [0.3 0.3 0.3 0.3]) ;
  R.Text1 = uicontrol(R.F,'style','text','units','normalized','string','Dimensions:','position',[0.1 0.80 0.25 0.085],'callback','','backgroundcolor',get(R.F,"Color")) ;
  R.Text2 = uicontrol(R.F,'style','text','units','normalized','string','Sorting   :','position',[0.1 0.60 0.30 0.085],'callback','','backgroundcolor',get(R.F,"Color")) ;
  R.Text3 = uicontrol(R.F,'style','text','units','normalized','string','Packing   :','position',[0.1 0.40 0.25 0.085],'callback','','backgroundcolor',get(R.F,"Color")) ;
  R.Dime = uicontrol(R.F,'style','edit','units','normalized','string','[150 , 150]','position',[0.37 0.80 0.25 0.085],'callback','') ;
  R.Sort = uicontrol(R.F,'style','popupmenu','units','normalized','string',{"Very Well","Well","Moderately Well","Moderately","Poorly","Very Poorly"},'position',[0.37 0.60 0.3 0.085],'callback',"") ;
  R.Pack = uicontrol(R.F,'style','edit','units','normalized','string','-0.1','position',[0.37 0.40 0.25 0.085],'callback','') ;
  R.Gene = uicontrol(R.F,'style','pushbutton','units','normalized','string','Generate','position',[0.35 0.1 0.35 0.085],'callback',@GetData) ;

  guidata (R.F, R) ;
  GetData(R.F,true)
end
