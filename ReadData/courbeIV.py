#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Ce module permet d'importer les données contenues dans un fichier I(V).
"""


#TODO : à supprimer



import sys
import os
try:
	import numpy
except ImportError:
    raise ImportError, 'Ce programme requiert la libraire Python-Scipy. http://www.scipy.org/'

def _line2floats(a):
   """_line2floats(string)
   Convertit une chaîne en liste de floats.
   Cette chaîne est typiquement une ligne d'un fichier de données :
   exemple : "1.25600E-1	2.3600E-3 ..."
   """
   return map(float,a.split())


def lireData(filename,format="array",col=4,paramcol=0,comment='V'):
	"""Vg,G = lireData(filename,format,col,comment)
	Retourne les données contenues dans un fichier de type I(V).
	Format du fichier d'entrée : colonnes Vg X Y T G
	
	arguments :
	    filename : nom du fichier à lire
	    format   : format des vecteurs retournés
	                  array : vecteur numpy (valeur par défaut)
	                  liste : liste python
	    col      : colonne des données
	    paramcol : colonne du paramètre variable (0 par défaut)
	"""
	f = open(filename, 'rb')
	data = f.readlines()
	f.close()
	data = filter(lambda x: x[0]!=comment,data)
	data = map(_line2floats,data)
	data = [elem for elem in data if elem != []]

	if format=="liste":
		Vg = [dat[paramcol] for dat in data]
		G = [dat[col] for dat in data]
	else :
		data=numpy.array(data)
		Vg=data[:,paramcol]
		G=data[:,col]
	return Vg,G

def selectData(Vg,Vgliste,*args):
	""" Vg2,G1,G2,... = select_donnees(Vg,Vgliste,G1,G2,...)
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
