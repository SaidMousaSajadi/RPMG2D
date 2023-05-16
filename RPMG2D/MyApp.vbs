' Dim fs, f
Set WshShell = CreateObject("WScript.Shell")
Set fs = CreateObject("Scripting.FileSystemObject")
Set f = fs.GetFile("MyApp.vbs")
Path = UCase(f.ParentFolder)
FileName = Path & "\MyApp.bat"
' MsgBox(FileName)
WshShell.Run chr(34) & FileName & Chr(34), 0
Set WshShell = Nothing
