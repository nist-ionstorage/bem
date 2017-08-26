import math

fileIn = "DieAssembly.xic"
fileOut = "print.eps"

size = 8*72 # Final size of postcript output page in pt
psdigits = 10.

###############################################
## Determine scaling factor from bounding box
###############################################
fxic = open(fileIn,"r")
firstloop = True
for line in fxic:
    args=line.rstrip().rstrip(';').split() #strip newline then split
    if args[0]=='P':
        # Extract coordinate pairs and scale
        coords = []
        for i in range((len(args)-1)/2):
            coords.append((int(args[2*i+1]),int(args[2*i+2])))
        # Find bounding box
        if firstloop:
            xboxmin = xboxmax = coords[0][0]
            yboxmin = yboxmax = coords[0][1]
            firstloop = False
        for x,y in coords:
            xboxmin = min(xboxmin,x)
            xboxmax = max(xboxmax,x)
            yboxmin = min(yboxmin,y)
            yboxmax = max(yboxmax,y)
fxic.close()
# Scale factor
scale = psdigits*size/float(max(xboxmax-xboxmin,yboxmax-yboxmin))

################################
## Determine shared edges
################################
fxic = open(fileIn,"r")
layer = 'none'
edgeList = {}
for line in fxic:
    args=line.rstrip().rstrip(';').split() #strip newline then split
    if args[0]=='L':
        layer = args[1]
    elif args[0]=='P':
        # Extract coordinate pairs and scale
        coords = []
        for i in range((len(args)-1)/2):
            coords.append((int(scale*(int(args[2*i+1])-xboxmin)),int(scale*(int(args[2*i+2])-yboxmin))))
        # Record all edges with their layer name.
        (xp,yp) = coords[len(coords)-1]
        for x,y in coords:
            key1 = (layer,xp,yp,x,y)
            key2 = (layer,x,y,xp,yp)
            if edgeList.has_key(key1):
                edgeList[key1] = 1   # Edge is shared
            elif edgeList.has_key(key2):
                edgeList[key2] = 1
            else:
                edgeList[key1] = 0   # Create key
            xp = x
            yp = y
            
fxic.close()


################################
## Output geometry to postscript
################################
fxic = open(fileIn,"r")
psfile = open(fileOut,"w")
psfile.write("%!\n")
psfile.write("%%%%BoundingBox: 0 0 %i %i\n" % (int((xboxmax-xboxmin)*scale/psdigits),int((yboxmax-yboxmin)*scale/psdigits)))
psfile.write("1 setlinewidth\n")
psfile.write("%f %f scale\n" % (1./psdigits, 1./psdigits))
layer = 'none'
for line in fxic:
    args=line.rstrip().rstrip(';').split() #strip newline then split
    if args[0]=='L':
        layer = args[1]
        if layer=='MTL1':
            color = ( 140./255., 10./255., 10./255. ) # RGB
            fillpattern = 1 # 0-none, 1-lines
            fillcolor = ( 0.8, 0.37, 0.1 ) # RGB 
            fillangle = -45
            fillspace = 0.05*psdigits*72
        elif layer=='MTL2':
            color = ( 30./255., 128./255., 46./255. ) # RGB
            fillpattern = 1 # 0-none, 1-lines
            fillcolor = ( 0.5, 0.8, 0 ) # RGB
            fillangle = 45
            fillspace = 0.05*psdigits*72
        elif layer=='OCUT':
            color = ( 46./255. , 30./255., 128./255. ) # RGB
            fillpattern = 1 # 0-none, 1-lines
            fillcolor = ( 0, 0.5, 0.8 ) # RGB
            fillangle = 45
            fillspace = 0.025*psdigits*72
        else:
            color = ( 0., 0., 0. ) # RGB
            fillpattern = 0 # 0-none, 1-lines
        psfile.write("%f %f %f setrgbcolor\n" % color)
    elif args[0]=='P':
        # Extract coordinate pairs and scale
        coords = []
        for i in range((len(args)-1)/2):
            coords.append((int(scale*(int(args[2*i+1])-xboxmin)),int(scale*(int(args[2*i+2])-yboxmin))))
        # Draw unshared boundary edges
        (xp,yp) = coords[len(coords)-1]
        for x,y in coords:
            key1 = (layer,xp,yp,x,y)
            key2 = (layer,x,y,xp,yp)
            if edgeList.has_key(key1) and edgeList[key1]==0:
                psfile.write("newpath %i %i moveto %i %i lineto stroke\n" \
                    % (xp,yp,x,y))
            elif edgeList.has_key(key2) and edgeList[key2]==0:
                psfile.write("newpath %i %i moveto %i %i lineto stroke\n" \
                    % (xp,yp,x,y))
            xp = x
            yp = y
        # Find bounding box for hatch fills
        xmin = xmax = coords[0][0]
        ymin = ymax = coords[0][1]
        for x,y in coords:
            xmin = min(xmin,x)
            xmax = max(xmax,x)
            ymin = min(ymin,y)
            ymax = max(ymax,y)
        # Draw hatching using a clipping path
        if fillpattern != 0:
            psfile.write("gsave %f %f %f setrgbcolor\n" % fillcolor)
            psfile.write("newpath %i %i moveto\n" % (coords[0][0], coords[0][1]))
            for x,y in coords[1:]:
                psfile.write("%i %i lineto\n" % (x, y))
            psfile.write("closepath clip\n")
        if fillpattern==1 and fillangle==-45:
            imin = int(math.floor((xmin+ymin)/fillspace/math.sqrt(2)))
            imax = int(math.ceil((xmax+ymax)/fillspace/math.sqrt(2)))
            for i in range(imin,imax+1):
                xa = int(i*fillspace*math.sqrt(2))-ymin
                xb = xa-(ymax-ymin)
                psfile.write("newpath %i %i moveto %i %i lineto stroke\n" \
                    % (xa,ymin,xb,ymax))
        elif fillpattern==1 and fillangle==45:
            imin = int(math.floor((xmin-ymax)/fillspace/math.sqrt(2)))
            imax = int(math.ceil((xmax-ymin)/fillspace/math.sqrt(2)))
            for i in range(imin,imax+1):
                xa = int(i*fillspace*math.sqrt(2))+ymax
                xb = xa-(ymax-ymin)
                psfile.write("newpath %i %i moveto %i %i lineto stroke\n" \
                    % (xa,ymax,xb,ymin))     
        if fillpattern != 0:
            psfile.write("grestore\n")
            
fxic.close()
psfile.close()
    