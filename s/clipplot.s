clip.plot <- function(width=10.5,height=8)
             dev.print(device=win.printer,format="metafile",
                       file="clipboard",width=width,height=height)
