#!/usr/bin/env python
"""
Generate a matplotlib colormap dictionary from a GMT cpt file.
Created on Sep 12, 2012

@author: behry
"""


def gmtColormap(fileName, GMTPath=None):
    import colorsys
    import numpy as N
    if type(GMTPath) == type(None):
        filePath = "/usr/local/cmaps/" + fileName + ".cpt"
    else:
        filePath = GMTPath + "/" + fileName + ".cpt"
    try:
        f = open(filePath)
    except:
        print "file ", filePath, "not found"
        return None

    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if l[0] == "#":
            if ls[-1] == "+HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    x = N.array(x , N.float)
    r = N.array(r , N.float)
    g = N.array(g , N.float)
    b = N.array(b , N.float)
    print colorModel
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr ; g[i] = gg ; b[i] = bb
    if colorModel == "RGB":
        r = r / 255.
        g = g / 255.
        b = b / 255.
    xNorm = (x - x[0]) / (x[-1] - x[0])
    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    return (colorDict)

if __name__ == '__main__':
    GMTPath = '/usr/local/GMT4.5.7/share/cpt/'
    filename = 'GMT_rainbow'
    cmap = gmtColormap(filename, GMTPath)
