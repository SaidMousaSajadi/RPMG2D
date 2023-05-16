function ModuleSelector() ;
  function GetData(obj , init = false)
    h = guidata(obj) ;
    setappdata(0,'GeneMethod','VCM') ;
    switch (gcbo)
      case {h.VCM}
        setappdata(0,'GeneMethod','VCM') ;
        close(h.F)
      case {h.VSM}
        setappdata(0,'GeneMethod','VSM') ;
        close(h.F)
      case {h.VMM}
        setappdata(0,'GeneMethod','VMM') ;
        close(h.F)
    endswitch
  end

  R.F = figure("toolbar", "none",'uicontextmenu',[],'menubar','none','name',"Method Selector To Generate Shapes",'NumberTitle','off','resize','off','units','normalized',"Position", [0.3 0.3 0.3 0.3]) ;
  R.Text1 = uicontrol(R.F,'style','text','units','normalized','string','Method Selector to Distribution Data of Random Shapes:','position',[0.05 0.80 0.85 0.085],'callback','','backgroundcolor',get(R.F,"Color")) ;
  R.VCM = uicontrol(R.F,'style','pushbutton','units','normalized','string','Variant Coefficient model(VCM) - Moderate Similar','position',[0.15 0.6 0.75 0.085],'callback',@GetData) ;
  R.VSM = uicontrol(R.F,'style','pushbutton','units','normalized','string','Variant Sigma model(VSM) - Very Similar','position',[0.15 0.4 0.75 0.085],'callback',@GetData) ;
  R.VMM = uicontrol(R.F,'style','pushbutton','units','normalized','string','Variant Mu model(VMM) - Poorly Similar','position',[0.15 0.2 0.75 0.085],'callback',@GetData) ;
  guidata (R.F, R) ;
  GetData(R.F,true)
end
