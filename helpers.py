# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.optimize import fsolve
from decimal import Decimal
import pandas as pd
from uncertainties import ufloat, unumpy
from uncertainties.umath import *
from subprocess import call, Popen, PIPE
import os

'''
    Fetch data from a csv file
    @param {string} file        - csv file containing data
    @param {string} col         - column header name
    @param {float|array} err    - error of the data set 
    @return {numpy.ndarray}		- numpy array containing the data (and their errors if given)
'''
def fetch (file, col, err=0):
	df = pd.read_csv(file)
	if err == 0:
		return np.array(df[col])
	else:
		return unumpy.uarray(df[col], np.zeros(df[col].size) + err)

def nominal (a):
	arr = np.array(a)
	for i, v in enumerate(arr):
		arr[i] = v.n if hasattr(v, 'n') else v
	return arr

def stddev (a):
	arr = np.array(a)
	for i, v in enumerate(arr):
		arr[i] = v.std_dev if hasattr(v, 'std_dev') else 0.0
	return arr

def equal (A, B, tol):
	for d in np.abs(A - B) - tol:
		if d > 0:
			return False
	return True

def compile (pdffile = "main", converter = "pdflatex"):
	#call([converter, '-jobname=' + pdffile, 'log/intermediate.tex'], stdin=None, stdout=None, stderr=None)
	p = Popen([converter, '-halt-on-error', '-jobname=' + pdffile, 'log/intermediate.tex'], stdin=None, stdout=PIPE, stderr=PIPE)
	output, err = p.communicate(b"input data that is passed to subprocess' None")
	rc = p.returncode

	print(rc)
	if rc != 0:
		print("Compilation of main.tex failed:")
		print(output)
	else:
		p = Popen(['evince', pdffile + '.pdf'], stdin=None, stdout=None, stderr=None)
		p.communicate(b"Opening pdf")
		os.remove(pdffile + ".aux")
		os.rename(pdffile + ".log", "log/" + pdffile + ".log")



'''
    Format a number (with its error if given) in latex
    @param {int|float|ufloat} nr	- number
    @param {string} sign			- format
    @return {string}				- formatted number in latex
'''
def fmt_number (nr, sign = '.2f'):

	s = "NaN"
	fmtNull = (r"{0:"+sign+r"}").format(0.0)
	if (np.abs(nr) > 9999 or np.abs(nr) < 0.1) and nr != 0:
		if hasattr(nr, 'n'):
			float_str = "{0:e}".format(nr.n)
			base, exponent = float_str.split("e")
			base, exponent = float(base), int(exponent)
			err = nr.s / 10**exponent
			if (nr.n == 10**exponent):
				s = (r"10^{{{0}}} \pm {1}").format(exponent, fmt_number(nr.s, sign))
			elif (r"{0:"+sign+r"}").format(err) == fmtNull:
				s = (r"{0:"+sign+r"} \times 10^{{{1}}} \pm {2}").format(base, exponent, fmt_number(nr.s, sign))
			else:
				s = (r"\left( {0:"+sign+r"} \pm {1:"+sign+r"} \right) \times 10^{{{2}}}").format(base, err, exponent)

		else:
			float_str = "{0:e}".format(nr)
			base, exponent = float_str.split("e")
			base, exponent = float(base), int(exponent)
			if (nr == 10**exponent):
				s = (r"10^{{{0}}}").format(exponent)
			else:
				s = (r"{0:"+sign+r"} \times 10^{{{1}}}").format(base, exponent)
	else:
		if hasattr(nr, 'n'):
			if (nr.s == 0.0):
				s = (r"{0:"+sign+r"}").format(nr.n)
			elif (r"{0:"+sign+r"}").format(nr.s) == fmtNull:
				s = (r"{0:"+sign+r"} \pm {1}").format(nr.n, fmt_number(nr.s, sign))
			else:
				s = (r"{0:"+sign+r"} \pm {1:"+sign+r"}").format(nr.n, nr.s)
		else:
			if isinstance(nr, int):
				s = str(nr)
			else:
				s = (r"{0:"+sign+r"}").format(nr)

	return s

def fmt_table (m, fmt=True):
	m = np.atleast_2d(m)
	t = []
	for i, v in enumerate(m):
		if fmt:
			t.append(' & '.join(format(x) for x in v))
		else:
			t.append(' & '.join(str(x) for x in v))
	
	tbl = ''
	tbl += t[0] + '\\\\\n\\hline\n\hline\n'
	t.pop(0)
	tbl += '\\\\\n'.join(str(x) for x in t)
	return tbl

def format(s, fmt=True):
	if isinstance(s, str): # it is a string
		return s
	elif hasattr(s, "__len__"): # it is an array
		return fmt_table(s, fmt)		
	elif (hasattr(s, 'n')): # it has uncertainties
		return "$" + fmt_number(s) + "$" if fmt else str(s)
	else: # it is a regular number (float, int)
		return "$" + fmt_number(s) + "$" if fmt else str(s)

def fitline(C, deg=None):
	if (deg == None):
		deg = len(C)-1
	pf = ""
	for i, c in (enumerate(C)):
		e = len(C) - (i+1) # exponent
		if e > deg:
			continue
		if c > 0 and deg != e:
			pf += "+"
		pf += fmt_number(c)
		if e == 1:
			pf += " x"
		elif e != 0:
			pf += " x^{" + str(len(C) - (i+1)) +"}"
	return pf


first = True
def replace(s, r, fmt=True, texfile = "main.tex"):
	global first
	if first:
		file = texfile
	else:
		file = 'log/intermediate.tex'
	# Read in the file
	with open(file, 'r') as file :
		filedata = file.read()

	# Replace the target string
	filedata = filedata.replace("[[" + s + "]]", format(r, fmt))

	# Write the file out again
	with open('log/intermediate.tex', 'w') as file:
		file.write(filedata)

	first = False

def upolyfit(x, y, deg, **kwargs):
	c, cov = np.polyfit(x, y, deg, cov=True, **kwargs)
	return unumpy.uarray(c, np.diagonal(cov))