#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Imports column-stored data from ascii files.
"""

import re
try:
    import numpy
except ImportError:
    raise ImportError, 'Ce programme requiert la libraire Python-Scipy. http://www.scipy.org/'


def float_or_zero(a):
   try:
      a_float = float(a)
   except ValueError:
      a_float = 0
   return a_float


def _line2floats(elems_raw):
   """_line2floats(string)
   Converts a string into a list of floats.
   The input string is typically like a line from a data file:
       example : "1.25600E-1   2.3600E-3 ..."
   """
   return map(float_or_zero, elems_raw.split())


def readFile(filename,cols=[],mode=None,format='array'):
    """X1,X2,... = readFile(filename,cols,comment,mode,format)
    Imports column-stored data from an ascii file.

    inputs :
        filename : file name (uh!)
        cols     : columns to be read (list of indexes).
        format   : format of returned data
                      'array' : numpy 1D-array (default value)
                      'list'  : python list
        mode     : if specified, define the list of columns to read
    """
    f = open(filename, 'rb')
    data = f.readlines()
    f.close()
    data = filter(lambda x: ((x[0]>='0') and (x[0]<='9')) or (x[0]=='-'), data)
    data = map(_line2floats,data)
    # remove empty lines:
    data = [elem for elem in data if elem != []]

    if mode is not None:
        if mode=='courbeIV':
            cols=[0,4]
        if mode=='rhombus':
            cols=[0,1,5]
        if mode=='rhombusDC':
            cols==[0,1,5,6]


    result = ()
    if format=='list':
        for col in cols:
            result += ( [dat[col] for dat in data] ,)
    else:
        data=numpy.array(data)
        for col in cols:
            result += ( data[:,col] ,)
    return result


def readDelftFile(filename, params, cols=[0,1]):
    f = open(filename, 'rb')
    
    data = []
    indices_2D_var = []
    vars={}
    sweep_var_value = 0
    re_line = re.compile("([\d|\-]+)\s+(.*)")
    re_var_line = re.compile("(\S+)=(\S+)(.*)")
    
    for line in f:
        a = re_line.match(line)
        if a is not None:
            lineID = int(a.group(1))
            line_data = a.group(2)
            if (lineID < 0) and (line_data[0] == 'L'):
                while True:
                    a = re_var_line.search(line_data)
                    if a is not None:
                        (var_name, var_value, line_data) = a.groups()
                        if not vars.has_key(var_name):
                            vars[var_name] = []
                        list_len = len(vars[var_name])
                        vars[var_name].append(var_value)
                        #vars[var_name] += [None for i in range(list_len+1+lineID)] + [var_value]
                    else: break
            if lineID > 0:
                data += [_line2floats(line_data)]
                indices_2D_var += [lineID-1]
    
    result = ()
    
    indices_2D_var = numpy.array(indices_2D_var)
    indices_2D_var -= indices_2D_var.min()
    sweep_var = params['sweep_var']
    if indices_2D_var[-1] > 0:
        if vars.has_key(sweep_var):
            result += ( numpy.array(map(float_or_zero, vars[sweep_var]))[indices_2D_var] ,)
        else:
            result += ( numpy.array(indices_2D_var) ,)
    
    """
    for (key, value) in vars.items():
    dico['col_var1'] = col_var1
    dico['col_var2'] = col_var2
    dico['var1_name'] = var1_name
    dico['var2_name'] = var2_name
    """
    
    if vars.has_key('sw'):
        params['ylabel'] = vars['sw'][-1]
    params['xlabel'] = sweep_var
    
    data=numpy.array(data)
    for col in cols:
        result += ( data[:,col] ,)
    return result


#TODO:réécrire cette fonction in a more general way :
#       vérifier que ca marche avec des listes ou des arrays
def selectData(Vg, Vgliste, *args):
    """ Vg2,G1,G2,... = select_donnees(Vg, Vgliste, G1, G2, ...)
    Retourne les données qui sont dans les fenêtres définies par la Vgliste
    [[Vgmin1,Vgmax1], ... , [VgminN, VgmaxN],...]
    """
    # on cherche le vecteur d'indice
    i = False
    for [a,b] in Vgliste:
        i = ( (Vg>=a) & (Vg<=b) ) | i
    # on fabrique une liste des vecteurs à découper
    result=[V]
    result.extend(args)
    return map(lambda w : w[i],result)


def coupe(V,v,*args,**kargs):
    """coupe(V,v,A,B,C,...,delta=1e-5) 
    Sélectionne les parties des vecteurs A, B, C,... qui correspondent
    à la valeur v (à +/- delta près) dans le vecteur V.
    exemples : coupe IdVg d'un rhombus pour Vd=vd 
               Vg,I = coupe(Vd,vd,Vg,I)

               coupe IdVd d'un rhombus pour Vg=vg 
               Vd,G,I = coupe(Vg,vg,Vd,G,I)
    """
    if kargs.has_key('delta'):
        delta = kargs['delta']
    else:
        delta = 1e-5
    i = numpy.nonzero(abs(V-v) <= delta)[0]
    if len(i)==0 :
        print "Attention, aucune donnée sélectionnée"
    return map(lambda w : w[i],args)
