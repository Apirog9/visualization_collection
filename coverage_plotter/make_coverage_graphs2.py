# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 15:44:18 2022

@author: Proccesing_PC2

make simple coverage graph for protein sequence, including quantified peptide count per position as well as rough fold change
quantitative difference p-value (uncorrected)
"""

from Bio import SeqIO
import re
from matplotlib import pyplot
import plotnine as plt
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

# FASTa file containing protein sequences
fastafile = "2021_07_13_uniprot-proteome_UP000005640.fasta"
# quantitation data file including peptide-level log transformed quantities per replicate and ID column compatible
# with ID in fasta file
precursordata = "combined_modified_peptide.tsv"
# names of columns with quantitation
conditions = ['pro_1','pro_2', 'pro_3', 'sen_1', 'sen_2', 'sen_3']
# names of columns of condition a for fold change calculation
cond_a = ['pro_1','pro_2', 'pro_3']
# names of columns of condition a for fold change calculation
cond_b = ['sen_1', 'sen_2', 'sen_3']
# minimum peptide number for protein to be plotted. Override names
minpeps = 1
# minimal number of valid values per condition for a peptide to be considered well-quantified
min_quants = 3
# file containing one identifier per line for proteins to be plotted
names = 'protgroups.list'
# column with protein name, should contain one identifier compatible with ID from FASTA file
idcolumn ='Protein ID'
# column with protein name, present in quantification dataframe, which will be used as plot title
namecolumn ='Entry Name'
# column with peptide sequences
peptidecol = 'Modified Sequence'
# use p or q value for color map use 'pval' or 'qval'
confidence = 'pval'
# if make log2 tranform before calculation for raw quantities
make_log2 = True






def makegraph_quant(sequence,pepdict,name,filename,confidence):  #make sure all peptides are used!!! not only quantitative
    labels = {'pval':'- log10 p-value','qval':'- log10 q-value'}
    cbar_label = labels[confidence]
    def getrange(minimum,maximum):
        # define clever FC axis range
        if abs(minimum - maximum) <2:
            min_new = minimum-1
            max_new = maximum+2
        return (min_new,max_new)
    # initiate point coordinates for line plot
    xdata = [x+1 for x in range(len(sequence))]
    ydata = [0 for x in range(len(sequence))]
    # initiate point coordinates for pseudo-line scatter plot
    xdata_quant = []
    ydata_quant = []
    pvals = []
    coverages = []
    # add quantitative data and coverage (detected peptide data) for respective lists
    for key in  pepdict.keys():
        peptide=key
        locs = [[x.start()+1,x.start()+len(peptide)] for x in re.finditer("=?"+peptide,sequence)]
        for loc  in locs:
            for point in range(loc[0],loc[1]):
                xdata_quant.append(point)
                ydata_quant.append(pepdict[key][0])
                pvals.append(pepdict[key][1])
        coverages = coverages+locs
    for coverage in coverages:
        ydata = ydata[:coverage[0]-1] + [x+1 for x in ydata[coverage[0]-1:coverage[1]]] + ydata[coverage[1]:]
        
        
    # create coverage (detected peptide count) line plot
    plot = plt.ggplot(plt.aes(xdata,ydata))
   
    plot = plot + \
    plt.labels.ggtitle(name) + \
    plt.scale_y_continuous(name = "AA detection count") + \
    plt.scale_x_continuous(name = "Position in sequence") + \
    plt.geom_line(plt.aes(x=xdata,y=ydata)) +\
    plt.theme(panel_background = plt.element_rect(fill ="white"),
              figure_size=(8, 4),
              plot_title=(plt.element_text(ha='center',size=18)),
              axis_title_x=(plt.element_text(ha='center',size=14,color='black')),
              axis_title_y=(plt.element_text(ha='center',size=14,color='black')),
              axis_text_y=(plt.element_text(ha='right',size=14,color='black')),
              axis_text_x=(plt.element_text(ha='center',size=14,color='black')))
    plot = plot.draw()
    # extract axes, add second axes and create pseudo-line scatter plot for quantitative values
    axes = plot.get_axes()
    pyplot.set_cmap("spring_r")
    axes2 =axes[0].twinx()
    im = axes2.scatter(xdata_quant,ydata_quant,marker = ".",alpha = 0.7, s=10,cmap='viridis_r', c = pvals,vmin=0,vmax=4)
    axes2.set_ylabel("Fold change",size=14)
    axes2.tick_params(axis='y',labelsize=14)

    if abs(min(ydata_quant) - max(ydata_quant)) <2:
        limits = getrange(min(ydata_quant),max(ydata_quant))
        axes2.set_ylim(limits)
    # modify colorbar paraneters
    cbar = pyplot.colorbar(im,pad=0.12,extend='neither')
    cbar.set_label(label=cbar_label,size=14)
    cbar.set_ticklabels(cbar.get_ticks(),size=14)
    cbar.solids.set(alpha=1)
    plot.savefig(filename+".png",dpi=1000)
    





"""
will set 1 to pval if NaN
will set 0 to foldchange in Exist_a and Exist_b False
will set foldchange to 2 if Valid_a only and Exist_b False will set pval to 0.05 
will set foldchange to -2 if Valid_b only and Exist_a False will set pval to 0.05 
"""
def filterrows_simple(serieslike):
    # assign False if all values in serieslike are NaN
    serieslike = list(serieslike)
    nans = len([x for x in serieslike if np.isnan(x)])
    
    if nans == len(serieslike):
        returnval = False
    else:
        returnval = True
    return returnval

def filterrows_minvals(serieslike,minvals):
    # assign False if minimum valid values in serieslike is not passed
    serieslike = list(serieslike)
    nans = len([x for x in serieslike if np.isnan(x)])
    
    if nans > len(serieslike) - minvals:
        returnval = False
    else:
        returnval = True
    return returnval

def welch(serieslike,a,b):
    # make simple Welch test, return p-value
    
    a_vals = list(serieslike[a])
    b_vals = list(serieslike[b])
    a_vals = [x for x in a_vals if not np.isnan(x)]
    b_vals = [x for x in b_vals if not np.isnan(x)]
    if len(a_vals)>2 and len(b_vals)>2:
        result = stats.ttest_ind(a_vals, b_vals,equal_var = False)
    else:
        result=1
    
    return result.pvalue

def foldchange_value(serieslike):
    returnval = None
    if serieslike["Exist_a"] == False or serieslike["Exist_b"] == False:
        returnval = 0
    if serieslike["Valid_a"] == True and serieslike["Exist_b"] == False:
        returnval = 2
    if serieslike["Exist_a"] == False and serieslike["Valid_b"] == True:
        returnval = -2
    if serieslike["Exist_a"] == True and serieslike["Exist_b"] == True:
        returnval = serieslike["a"] - serieslike["b"]
    
    return returnval

def welch_coloring(serieslike,welch_q_val):
    # returns approximate and arbitrary values enabling visualization of incompletely quantified peptides
    # peptides with one condition valid and one empty gets q=0.01
    # peptide with no valid quantification gets q=1
    # if q val is wrongly asigned as 0 (e.g using externally done test) recieves q= 0.0001
    #if q val is wrongly asigned as NaN (e.g using externally done test) recieves q= 1
    
    
    returnval = None
    if serieslike["Exist_a"] == False or serieslike["Exist_b"] == False:
        returnval = 1
    if serieslike["Valid_a"] == True and serieslike["Exist_b"] == False:
        returnval = 0.01
    if serieslike["Exist_a"] == False and serieslike["Valid_b"] == True:
        returnval = 0.01
    if serieslike["Exist_a"] == True and serieslike["Exist_b"] == True:
        if np.isnan(serieslike[welch_q_val]):
            returnval = 1
        elif serieslike[welch_q_val] == 0:
                returnval = 0.0001
        else:
            returnval = serieslike[welch_q_val]

    returnval = -1* np.log10(returnval)

    return returnval
    

def graph_func_quant(protname,namecolumn,peptidecolumn,dataframe,conditions,fastadict,a,b,confidence):

    
    if dataframe.shape[0] >= minpeps:
        name = dataframe[namecolumn].iloc[0]
        for item in fastadict.keys():
            if protname ==item.split("|")[1]:
                seq = str(fastadict[item].seq)
                pepsdict = {}
                for row in dataframe.iterrows():

                    pepsdict[str(row[1][peptidecolumn])] = [float(row[1]["Fold_change"]),float(row[1]["Welch_p_value_coloring"])]

                result = makegraph_quant(seq,pepsdict,name,protname,confidence)
  
    else:
        pass
    return result

def quantitate_rough(dataframe,conditions,a,b,minvals,confidence='pval'):
    coldict = {'pval':'Welch_p_value','qval':'Welch_p_value'}
    color_column = coldict[confidence]
    dataframe["Exist"] = dataframe[conditions].apply(filterrows_simple,axis=1)
    dataframe = dataframe[dataframe["Exist"] == True]
    dataframe["Exist_a"] = dataframe[a].apply(filterrows_simple,axis=1)
    dataframe["Exist_b"] = dataframe[b].apply(filterrows_simple,axis=1)
    dataframe["Valid_a"] = dataframe[a].apply(filterrows_minvals,axis=1,args = [minvals])
    dataframe["Valid_b"] = dataframe[b].apply(filterrows_minvals,axis=1,args = [minvals])
    dataframe["a"] = dataframe[cond_a].mean(skipna=True,axis=1)
    dataframe["b"] = dataframe[cond_b].mean(skipna=True,axis=1)
    dataframe["Fold_change"] = dataframe[["a","b","Exist_a","Exist_b","Valid_a","Valid_b"]].apply(foldchange_value,axis=1)
    #dataframe["Welch_p_value"] = np.nan
    dataframe["Welch_p_value"] = dataframe[conditions].apply(welch,axis=1,args = [a,b])
    print(list(dataframe["Welch_p_value"]))
    dataframe['Welch_q_value'] = multipletests(dataframe["Welch_p_value"],method = 'fdr_bh')[1]
    dataframe["Welch_p_value_coloring"] = dataframe.apply(welch_coloring,axis=1,args = [color_column])
    
    return dataframe







precursordata = pd.read_csv(precursordata, sep = "\t")
print(precursordata.columns)

if make_log2:
    precursordata[conditions] = np.log2(precursordata[conditions])
precursordata = quantitate_rough(precursordata,conditions,cond_a,cond_b,min_quants,confidence)
names = open(names,'r').read().splitlines()
precursordata_group = precursordata.groupby(idcolumn)
fastadict = SeqIO.to_dict((SeqIO.parse(fastafile, "fasta")))


for name in names:
    try:
        tempframe = precursordata_group.get_group(name)
        tempframe = pd.DataFrame(tempframe)
        a = graph_func_quant(name,namecolumn,peptidecol,tempframe,conditions,fastadict,cond_a,cond_b,confidence)
    except KeyError:
        print('no identifier '+name)
    
    
    
    






