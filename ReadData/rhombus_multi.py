#!/usr/bin/python
# -*- coding: utf-8 -*-

#
"""
Import and plot in color scale 2D data.
"""

import StringIO
import Image
import warnings
import matplotlib.pyplot
import pyx

# intra-package modules
import readFile
import palette
import constants

try:
    import numpy
except ImportError:
    raise ImportError, 'Ce programme requiert la libraire Python-numpy. http://www.numpy.org/'

#pyx.text.set(mode='latex')

class Data:
    """Store 2D data"""
    def __init__(self,data,Nx,Ny,Xmin,Xmax,Ymin,Ymax,gmin,gmax):
        self.data=data
        self.Nx=Nx
        self.Ny=Ny
        self.Xmin=Xmin
        self.Xmax=Xmax
        self.Ymin=Ymin
        self.Ymax=Ymax
        self.gmin=gmin
        self.gmax=gmax

def loadRhombus(filename, measurescale=1.0/constants.G0_SI, paramscale=[1,1], col=5, paramcol=[0,1], stepX=None, stepY=None):
    """data=loadRhombus(filename)
    Load column-stored 2D data as a matrix. Each line of the file corresponds to a data point.
    Output: an instance of class Data.
    Optional inputs:
        col          : column index of the data to read.
        paramcol     : list of two index corresponding to the two parameters [X,Y]
        measurescale : the raw data are multiplicated by this factor
        paramscale   : the parameters are multiplicated by these factors [scaleX,scaleY]
        stepX,stepY  : step of variation of the parameters.
    """
    Vg,Vd,G = readFile.readFile(filename,cols=paramcol+[col])
    G = G * measurescale
    Vg = Vg*paramscale[0]
    Vd = Vd*paramscale[1]
    return makeMatrix(Vg,Vd,G,stepX=stepX,stepY=stepY)

def loadRhombusmulti(filename, Yfirst, stepY, N, Nfirst=0, measurescale=1.0/constants.G0_SI, paramscale=1, cols=[0,4], stepX=None):
    """data=loadRhombus(filename)
    Load column-stored 2D data as a matrix. Each line of the file corresponds to a data point.
    Output: an instance of class Data.
    Optional inputs:
        col          : column index of the data to read.
        paramcol     : list of two index corresponding to the two parameters [X,Y]
        measurescale : the raw data are multiplicated by this factor
        paramscale   : the parameters are multiplicated by these factors [scaleX,scaleY]
        stepX,stepY  : step of variation of the parameters.
    """
    Vg=numpy.zeros((0,1),dtype=float)
    Vd=numpy.zeros((0,1),dtype=float)
    G=numpy.zeros((0,1),dtype=float)
    vd=Yfirst
    for i in range(Nfirst,Nfirst+N):
        vg,g = readFile.readFile(filename % i,cols=cols)
        Vg=numpy.concatenate((Vg,vg[:,numpy.newaxis]))
        G=numpy.concatenate((G,g[:,numpy.newaxis]))
        Vd=numpy.concatenate((Vd,vd*numpy.ones((len(vg),1))))
        vd = vd + stepY
    G = G * measurescale
    vg = Vg*paramscale
    return makeMatrix(Vg,Vd,G,stepX=stepX,stepY=abs(stepY))


def makeMatrix(X,Y,data,stepX=None,stepY=None):
    """Data = makeMatrix(X,Y,data)
    X, Y and data are 1D-arrays.
    The data are rearranged into a 2D-array.
    An instance of the class Data is returned.
    """
    # On cherche le nombre de points dans X et Y
    x = numpy.unique(X)
    y = numpy.unique(Y)
    if stepX is None:
        stepX = x[1]-x[0]
        if stepX < 1e-9:
            warnings.warn("The X step seems very low. You should consider setting it by hand.") 
    if stepY is None:
        stepY = y[1]-y[0]
        if stepY < 1e-9:
            warnings.warn("The Y step seems very low. You should consider setting it by hand.")
    Nx = int(round((x.max()-x.min())/stepX,0))+1
    Ny = int(round((y.max()-y.min())/stepY,0))+1

    # On convertit les coordonnées de chaque point en indice dans la matrice finale    
    Ix = numpy.floor((X-x.min())/stepX+0.1).astype(int)
    Iy = numpy.floor((Y-y.min())/stepY+0.1).astype(int)

    # On écrit les données dans une matrice
    matrix = numpy.zeros((Ny,Nx))
    matrix[(Iy,Ix)] = data

    return Data(matrix,Nx,Ny,x.min(),x.max(),y.min(),y.max(),data.min(),data.max())


def dataselectrange(d,xmin=None,xmax=None,ymin=None,ymax=None):
    """select a sub-matrix in data"""
    if xmin is None:
        xmin=d.Xmin
    if xmax is None:
        xmax=d.Xmax
    if ymin is None:
        ymin=d.Ymin
    if ymax is None:
        ymax=d.Ymax
    Ixmin= max(int(round((xmin-d.Xmin)/(d.Xmax-d.Xmin)*(d.Nx-1))),0)
    Ixmax= min(int(round((xmax-d.Xmin)/(d.Xmax-d.Xmin)*(d.Nx-1))),d.Nx-1)
    Iymin= max(int(round((ymin-d.Ymin)/(d.Ymax-d.Ymin)*(d.Ny-1))),0)
    Iymax= min(int(round((ymax-d.Ymin)/(d.Ymax-d.Ymin)*(d.Ny-1))),d.Ny-1)
    xmin = d.Xmin + (d.Xmax-d.Xmin)/(d.Nx-1)*Ixmin
    xmax = d.Xmin + (d.Xmax-d.Xmin)/(d.Nx-1)*Ixmax
    ymin = d.Ymin + (d.Ymax-d.Ymin)/(d.Ny-1)*Iymin
    ymax = d.Ymin + (d.Ymax-d.Ymin)/(d.Ny-1)*Iymax
    #TODO : recalculer gmin et gmax (impossible à faire exactement puisque on est déjà sous forme matrice)
    return Data(d.data[Iymin:Iymax+1,Ixmin:Ixmax+1],Ixmax-Ixmin+1,Iymax-Iymin+1,xmin,xmax,ymin,ymax,d.gmin,d.gmax)

def datatranslate(d,x=0.0,y=0.0):
    """translate the data into the X-Y plan"""
    return Data(d.data,d.Nx,d.Ny,d.Xmin+x,d.Xmax+x,d.Ymin+y,d.Ymax+y,d.gmin,d.gmax)

def datatranspose(d):
    return Data(d.data.transpose(),d.Ny,d.Nx,d.Ymin,d.Ymax,d.Xmin,d.Xmax,d.gmin,d.gmax)

def normalise(data,gmin=None,gmax=None,zero=0.0,top=1.0,down=0.0,mid=0.5):
    """G = normalise(data)
    normalise the data:
        gmax --> top
        zero --> mid    (omitted if zero is set to None)
        gmin --> down
    default values:
        gmax max(data)
        zero 0.0
        gmin min(data)
        top  1.0
        mid  0.5
        down 0.0
    """
    G=numpy.zeros(data.shape)
    if gmin is None:
        gmin=data.min()
    if gmax is None:
        gmax=data.max()
    if zero is not None:
        ip = numpy.nonzero(data>zero)
        im = numpy.nonzero(data<zero)
        iz = numpy.nonzero(data==zero)
        G[ip]= mid+(top-mid)*(data[ip]-zero)/(gmax-zero)
        G[im]= mid+(down-mid)*(data[im]-zero)/(gmin-zero)
        G[iz]= mid
    else:
        G=down + (top-down)*(data-gmin)/(gmax-gmin)
    return G

def logscale(data, seuil=1.):
    """G=logscale(data,seuil)
    apply sgn(.)*log((abs(./seuil)).
    Do nothing is seuil is None or 0.
    """
    if (seuil is None) or (seuil is 0):

        G = data
    else:
        mask = abs(data) <= seuil
        G = numpy.sign(data)*numpy.log10(abs(data)/seuil)
        G[mask] = 0
        del mask
    
    return G
def logscale_reverse(data, seuil=1.):
    return seuil*10**data

def plot(d,gmin=None,gmax=None,zero=0,glog=None,palette=palette.rwb,width=8,height=6,space = 0.05,colorbarwidth = 0.05,xlabel=r'$V_g$ (mV)',ylabel=r'$V_d$ (mV)',glabel=r'G $(\frac{e^2}{h})$',insert_graph=True,insert_colorbar=True,rotate=True,top=1.0,down=0.0,mid=0.5,
xmin=None,xmax=None,ymin=None,ymax=None,title=None):
    """c,g=plot(data)
        make a 2D-colorplot of data (a Data class instance)
        return c : a python pyx canvas instance
               g : a python pyx graph instance
    optionnal parameters :
    - data scaling :
        glog : parameter for logarithmic scaling of data
        gmax --> top
        gmin --> down
        zero --> mid  (if zero is not None)
        palette
    - xmin,xmax,ymin,ymax
    - formatting :
        width
        height
        space : space between graph and colorbar (% of width)
        colorbarwidth (% of width)
        xlabel, ylabel, glabel, title
        insert_graph : insert or not the graph g in canvas c (useful to work further on graph g)
        insert_colorbar
        rotate : set glabel orientation
    """
    if isinstance(d,tuple)==False:
        d=(d,)

    # looking for the maxima and minima in data, axis, ...
    if gmax is None:
        gmax=max( [di.data.max() for di in d])
    if gmin is None:
        gmin=min( [di.data.min() for di in d])

    # size of graph
    if xmin is None:
        xmin=min( [di.Xmin for di in d])
    if xmax is None:
        xmax=max( [di.Xmax for di in d])
    if ymin is None:
        ymin=min( [di.Ymin for di in d])
    if ymax is None:
        ymax=max( [di.Ymax for di in d])

    stepX=min( [(di.Xmax-di.Xmin)/(di.Nx-1.0) for di in d])
    Nx = (xmax-xmin)/stepX + 1
    stepY=min( [(di.Ymax-di.Ymin)/(di.Ny-1.0) for di in d])
    Ny = (ymax-ymin)/stepY + 1

    # to improve the color interpolation when exporting into PYX
    dx= (xmax-xmin-1)/Nx
    dy= (ymax-ymin-1)/Ny

    # figure
    dpi=300
    fig = matplotlib.pyplot.figure(figsize=(Nx*1.0/dpi, Ny*1.0/dpi),dpi=dpi, frameon=False)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    for di in d:
        if glog is None:
            cmin=gmin
            cmax=gmax
            G=di.data
        else:
            G = logscale(di.data,glog)
            if zero is not None:
                [cmax,cmin,zero] = logscale(numpy.array([gmax,gmin,zero]),glog)
            else:
                [cmax,cmin] = logscale(numpy.array([gmax,gmin]),glog)
        G = normalise(G,gmin=cmin,gmax=cmax,zero=zero,top=top,down=down,mid=mid)
        ax.imshow(G, extent=(di.Xmin,di.Xmax,di.Ymin,di.Ymax), cmap=palette, vmin=0,vmax=1, aspect='auto', origin='lower',interpolation='nearest')
    ax.set_xlim(xmin,xmax+dx)
    ax.set_ylim(ymin,ymax+dy)
    figfile = StringIO.StringIO()
    fig.savefig(figfile,format='png',dpi=dpi)
    matplotlib.pyplot.close(fig)

    # colorbar
    dat=numpy.linspace(gmin,gmax,250).reshape((250,1))
    if glog is not None:
        dat = logscale(dat,glog)
    dat = normalise(dat,gmin=cmin,gmax=cmax,zero=zero,top=top,down=down,mid=mid)
    colorbar = matplotlib.pyplot.figure(figsize=(10.0/dpi, 250.0/dpi),dpi=dpi, frameon=False)
    ax = colorbar.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.imshow(dat,cmap=palette,vmin=0,vmax=1,aspect='auto',origin='lower')
    colorbarfile = StringIO.StringIO()
    colorbar.savefig(colorbarfile,format='png',dpi=dpi)
    matplotlib.pyplot.close(colorbar)

    # rewind the data
    figfile.seek(0)
    colorbarfile.seek(0)
    fig = Image.open(figfile)
    fig = fig.convert('RGB')
    colorbar = Image.open(colorbarfile)
    colorbar = colorbar.convert('RGB')
    figfile.close()
    colorbarfile.close()

    c = pyx.canvas.canvas()
    c.insert(pyx.bitmap.bitmap(0, 0, fig, compressmode="Flate"),[pyx.trafo.scale(sx=width/pyx.unit.tocm(Nx*1.0/dpi*pyx.unit.t_inch),sy=height/pyx.unit.tocm(Ny*1.0/dpi*pyx.unit.t_inch))])

    g = pyx.graph.graphxy(width=width,height=height,
                x=pyx.graph.axis.linear(min=xmin, max=xmax, title=xlabel),
                y=pyx.graph.axis.linear(min=ymin, max=ymax, title=ylabel))

    #colorbar
    if insert_colorbar==True:
        c.insert(pyx.bitmap.bitmap(width*(1+space), 0, colorbar, compressmode="Flate"),[pyx.trafo.scale(sx=width*colorbarwidth/pyx.unit.tocm(10.0/dpi*pyx.unit.t_inch),sy=height/pyx.unit.tocm(250.0/dpi*pyx.unit.t_inch),x=width*(1+space),y=0)])
        if rotate==True:
            painter = pyx.graph.axis.painter.regular(titledirection=pyx.graph.axis.painter.rotatetext(270))
        else:
            painter = pyx.graph.axis.painter.regular(titledirection=pyx.graph.axis.painter.rotatetext(90))
        h = pyx.graph.graphxy(width=width*colorbarwidth,height=height,xpos=width*(1+space),
                    x=pyx.graph.axis.linear(min=0, max=1, parter=None),
                    y2=pyx.graph.axis.linear(min=gmin, max=gmax, title=glabel, painter=painter))
        c.insert(h)

    if title is not None:
        c.text(g.xpos+0.5*width, g.ypos+height*(1+space), title, [pyx.text.halign.boxcenter, pyx.text.valign.bottom])

    if insert_graph==True:
        c.insert(g)
    return c,g
