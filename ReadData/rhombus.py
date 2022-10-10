#!/usr/bin/python
# -*- coding: utf-8 -*-

#TODO refaire l'intro
"""
Ce module permet d'importer les données contenues dans un fichier rhombus, et d'extraire des coupes I(Vg) ou I(Vd).
"""
from matplotlib import rcParams
import StringIO
import Image
import warnings
import matplotlib.pyplot
import re
#import pyx

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
    def __init__(self,data,Nx,Ny,Xmin,Xmax,Ymin,Ymax):
        self.data=data
        self.Nx=Nx
        self.Ny=Ny
        self.Xmin=Xmin
        self.Xmax=Xmax
        self.Ymin=Ymin
        self.Ymax=Ymax

def loadRhombus(filename, params, measurescale=1.0/constants.G0_SI, paramscale=[1,1], col=5, paramcol=[0,1], stepX=None, stepY=None, Delft=False):
    if Delft is False:
        Vg,Vd,G = readFile.readFile(filename, cols=paramcol+[col])
    else:
        Vg,Vd,G = readFile.readDelftFile(filename, params, cols=[0, col])
    print "raw data: gmin =", G.min(), "; gmax =", G.max()
    G = G * measurescale - params['offset']
    Vg = Vg*paramscale[0]
    Vd = Vd*paramscale[1]
    return makeMatrix(Vg, Vd, G, stepX=stepX, stepY=stepY)


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
    Nx = int(round((x[-1]-x[0])/stepX,0))+1
    Ny = int(round((y[-1]-y[0])/stepY,0))+1

    # On convertit les coordonnées de chaque point en indice dans la matrice finale    
    Ix = numpy.around((X-x[0])/stepX).astype(int)
    Iy = numpy.around((Y-y[0])/stepY).astype(int)

    # On écrit les données dans une matrice
    matrix = numpy.zeros((Ny,Nx))
    matrix[(Iy,Ix)] = data

    return Data(matrix,Nx,Ny,x.min(),x.max(),y.min(),y.max())


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
    return Data(d.data[Iymin:Iymax+1,Ixmin:Ixmax+1],Ixmax-Ixmin+1,Iymax-Iymin+1,xmin,xmax,ymin,ymax)

def datatranslate(d,x=0.0,y=0.0):
    return Data(d.data,d.Nx,d.Ny,d.Xmin+x,d.Xmax+x,d.Ymin+y,d.Ymax+y)

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
        ip = data > zero
        im = data < zero
        iz = data == zero
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
    return numpy.sign(data)*seuil*10**numpy.abs(data)

def plot(d, gmin=None, gmax=None, zero=0, glog=None, palette=palette.rwb, width=8, height=6, space = 0.05, colorbarwidth = 0.05,
        xlabel=r'$V_g$ (mV)', ylabel=r'$V_d$ (mV)', glabel=r'G $(\frac{e^2}{h})$', insert_graph=True, insert_colorbar=True,
        rotate=True, top=1.0, down=0.0, mid=0.5, xmin=None, xmax=None, ymin=None, ymax=None):
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
        xlabel, ylabel, glabel
        insert_graph : insert or not the graph g in canvas c (useful to work further on graph g)
        insert_colorbar
        rotate : set glabel orientation
    """
    if gmax is None:
        gmax=d.data.max()
    if gmin is None:
        gmin=d.data.min()

    if glog is None:
        cmin=gmin
        cmax=gmax
        G=d.data
    else:
        G = logscale(d.data,glog)
        if zero is not None:
            [cmax,cmin,zero] = logscale(numpy.array([gmax,gmin,zero]),glog)
        else:
            [cmax,cmin] = logscale(numpy.array([gmax,gmin]),glog)
    print gmin
    G = normalise(G,gmin=cmin,gmax=cmax,zero=zero,top=top,down=down,mid=mid)

    # size of graph
    if xmin is None:
        xmin = d.Xmin
    if xmax is None:
        xmax = d.Xmax
    if ymin is None:
        ymin = d.Ymin
    if ymax is None:
        ymax = d.Ymax

    Nx = d.Nx * (xmax-xmin)/(d.Xmax-d.Xmin)
    Ny = d.Ny * (ymax-ymin)/(d.Ymax-d.Ymin)
    dx = (xmax-xmin-1)/Nx
    dy = (ymax-ymin-1)/Ny
    
    # figure
    dpi=1.
    fig = matplotlib.pyplot.figure(figsize=(Nx*1.0/dpi, Ny*1.0/dpi), dpi=dpi, frameon=False)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    matplotlib.pyplot.imshow(G, cmap=palette, extent=(d.Xmin,d.Xmax,d.Ymin,d.Ymax), aspect='auto', origin='lower', vmin=0.0, vmax=1.0)
    matplotlib.pyplot.xlim(xmin, xmax+dx)
    matplotlib.pyplot.ylim(ymin, ymax+dy)
    figfile = StringIO.StringIO()
    matplotlib.pyplot.savefig(figfile, format='png', dpi=dpi)
    #matplotlib.pyplot.savefig("figure.png", format='png', dpi=dpi)
    matplotlib.pyplot.close(fig)

    # colorbar
    dat=numpy.linspace(gmin,gmax,250).reshape((250,1))
    if glog is not None:
        dat = logscale(dat,glog)
    dat = normalise(dat,gmin=cmin,gmax=cmax,zero=zero,top=top,down=down,mid=mid)
    colorbar = matplotlib.pyplot.figure(figsize=(10.0/dpi, 250.0/dpi), dpi=dpi, frameon=False)
    ax = colorbar.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    matplotlib.pyplot.imshow(dat,cmap=palette,vmin=0,vmax=1,aspect='auto',origin='lower')
    colorbarfile = StringIO.StringIO()
    matplotlib.pyplot.savefig(colorbarfile,format='png',dpi=dpi)
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
    if insert_graph==True:
        c.insert(g)
    return c,g


def get_default_config():
    dico={}
    
    class ARG:
        def __init__(self, variable, type, default):
            self.variable = variable
            self.type = type
            self.default = default
    
    args=[]
    
    # Floats
    float_args = ['glog', 'gmin', 'gmax', 'xmin', 'xmax', 'ymin', 'ymax', 'ratioxy', 'stepX', 'stepY', 'col', 
        'font_size', 'xtick_major_pad', 'ytick_major_pad', 'dpi',
        'Mxx', 'Mxy', 'Myx', 'Myy', 'Vd',
        'triangles_point1_x', 'triangles_point1_y', 'triangles_point2_x', 'triangles_point2_y',
        'triangles_point3_x', 'triangles_point3_y', 'triangles_vect_x', 'triangles_vect_y']
    for arg in float_args:
        args.append(ARG(arg, 'float', None))
    args.append(ARG('paramscale_x', 'float', 1))
    args.append(ARG('paramscale_y', 'float', 1))
    args.append(ARG('measurescale', 'float', 1))
    args.append(ARG('offset', 'float', 0))
    
    # Strings
    string_args = ['data_label', 'colorbar_label', 'font']
    for arg in string_args:
        args.append(ARG(arg, 'string', ""))
    args.append(ARG('xlabel', 'string', None))
    args.append(ARG('ylabel', 'string', None))
    args.append(ARG('sweep_var', 'string', None))
    args.append(ARG('out_filename', 'string', None))
    args.append(ARG('mode', 'string', '2'))
    args.append(ARG('pref_quantity', 'string', 'cond'))
    
    # Bools
    boolean_args = ['gsym', 'neg_conductance', 'colorbar_display', 'usetex', 'Deflt_file']
    for arg in boolean_args:
        args.append(ARG(arg, 'bool', False))
    
    # Float lists
    float_list_args = ['ax_yticks', 'ax_xticks', 'colorbar_ticks']
    for arg in float_list_args:
        args.append(ARG(arg, 'float_list', None))
    
    # String lists
    string_list_args = ['colorbar_ticklabels']
    for arg in string_list_args:
        args.append(ARG(arg, 'string_list', None))
    
    for arg in args:
        dico[arg.variable] = arg.default
    
    return dico


def plot_bis(d, params, zero=0, palette=palette.rwb, width=8, height=6, space = 0.05, colorbarwidth = 0.05, 
        glabel=r'G $(\frac{e^2}{h})$', insert_graph=True, insert_colorbar=True,
        rotate=True, top=1.0, down=0.0, mid=0.5, fig=None, display_frame=True):
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
        xlabel, ylabel, glabel
        insert_graph : insert or not the graph g in canvas c (useful to work further on graph g)
        insert_colorbar
        rotate : set glabel orientation
    """
    
    matplotlib.rc('text', usetex=params['usetex'])
    if params['font_size'] is not None:
        matplotlib.rc('font', size=params['font_size'])
    if params['xtick_major_pad'] is not None:
        matplotlib.rcParams['xtick.major.pad'] = params['xtick_major_pad']
    if params['ytick_major_pad'] is not None:
        matplotlib.rcParams['ytick.major.pad'] = params['ytick_major_pad']
    if params['font'] != "":
        rcParams['font.family']='serif'
        rcParams['font.serif']=[params['font']]
    gmin=params['gmin']
    gmax=params['gmax']
    glog=params['glog']
    xmin=params['xmin']
    xmax=params['xmax']
    ymin=params['ymin']
    ymax=params['ymax']
    
    ax_xticks = params['ax_xticks']
    ax_yticks = params['ax_yticks']
    xlabel = params['xlabel']
    ylabel = params['ylabel']
    
    colorbar_display = params['colorbar_display']
    colorbar_label = params['colorbar_label']
    colorbar_ticks = params['colorbar_ticks']
    colorbar_ticklabels = params['colorbar_ticklabels']
    
    if gmax is None:
        gmax=d.data.max()
    if gmin is None:
        gmin=d.data.min()
    if params['gsym'] == True:
        gmax = max(abs(gmax), abs(gmin))
        gmin = -1*max(abs(gmax), abs(gmin))
    
    # colorbar
    cbar_values = numpy.linspace(gmin, gmax, 99)
    
    if glog is None:
        cmin=gmin
        cmax=gmax
        G=d.data
    else:
        G = logscale(d.data, glog)
        cbar_values = logscale(cbar_values, glog)
        if zero is not None:
            [cmax,cmin,zero] = logscale(numpy.array([gmax,gmin,zero]), glog)
        else:
            [cmax,cmin] = logscale(numpy.array([gmax,gmin]), glog)
    print "plot: gmin =", gmin, "; gmax =", gmax
    G = normalise(G, gmin=cmin, gmax=cmax, zero=zero, top=top, down=down, mid=mid)
    
    # colorbar
    cbar_values = normalise(cbar_values, gmin=cmin, gmax=cmax, zero=zero, top=top, down=down, mid=mid)
    cbar_Npoints_n = numpy.floor(100*(mid-down)/(top-down))
    if cbar_Npoints_n < 1:
        cbar_Npoints_n = 1
    if zero is None:
        cbar_boundaries = numpy.linspace(cmin, cmax, 99)
        if glog is not None:
            cbar_boundaries = logscale_reverse(cbar_boundaries, glog)
    else:
        cbar_boundaries_n = numpy.linspace(cmin, zero, cbar_Npoints_n)[:-1]
        cbar_boundaries_p = numpy.linspace(zero, cmax, 100-cbar_Npoints_n)
        if glog is not None:
            cbar_boundaries_n = logscale_reverse(cbar_boundaries_n, glog)
            cbar_boundaries_p = logscale_reverse(cbar_boundaries_p, glog)
        cbar_boundaries = numpy.append(cbar_boundaries_n, cbar_boundaries_p)
    
    # size of graph
    if xmin is None:
        xmin = d.Xmin
    if xmax is None:
        xmax = d.Xmax
    if ymin is None:
        ymin = d.Ymin
    if ymax is None:
        ymax = d.Ymax
    
    # figure
    if params['dpi'] is not None:
        dpi=params['dpi']
    else: dpi=72.
    
    if fig is None:
        if display_frame is False:
            fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=dpi, frameon=False)
            ax = fig.add_axes([0.0, 0.0, 1.0, 0.99], frameon=False)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=dpi, frameon=True)
            if colorbar_display:
                axes = [0.15, 0.2, 0.7, 0.7]
                matplotlib.pyplot.axes(axes)
            if params['axes'] is not None:
                axes = params['axes']
                matplotlib.pyplot.axes(axes)
    matplotlib.pyplot.imshow(gmin+G*(gmax-gmin), cmap=palette, extent=(d.Xmin,d.Xmax,d.Ymin,d.Ymax), aspect='auto', origin='lower', vmin=gmin, vmax=gmax, interpolation="nearest")
    matplotlib.pyplot.xlim(xmin, xmax)
    matplotlib.pyplot.ylim(ymin, ymax)
    matplotlib.pyplot.axis((xmin, xmax, ymin, ymax))
    matplotlib.pyplot.xlabel(xlabel, labelpad=10)
    matplotlib.pyplot.ylabel(ylabel, labelpad=10)
    if ax_xticks is not None:
        matplotlib.pyplot.xticks(ax_xticks)
    if ax_yticks is not None:
        matplotlib.pyplot.yticks(ax_yticks)
    
    # colorbar
    if colorbar_display:
        cb = matplotlib.pyplot.colorbar(boundaries=cbar_boundaries, values=numpy.linspace(gmin,gmax,98), format=None)
        cb.set_ticks(colorbar_ticks)
        cb.update_ticks()
        if colorbar_ticklabels is not None:
            cb.set_ticklabels(colorbar_ticklabels)
        if colorbar_label is not None:
            cb.set_label(colorbar_label, labelpad=10)
    
    return fig
