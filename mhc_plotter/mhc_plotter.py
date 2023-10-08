# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 20:41:22 2020

@author: APirog

Parsing and plotting NetMHCpan text result.
Peptidereuslt lists are also useful for annotation and storing data in more convinient formats

"""

import pandas as pd
import plotnine as plt


# here specify netmhcpan output files
netmhcpans = ["example_netmhcpan.txt"]



# use 0 to plot only binders. Total peptide number is user-supplied, because rarely all of the detected population 
# was used for NetMHCpan prediction. Alternatively, total peptide number from prediction is printed in console
total_quantified =  2406
# maximum percentrank for weak binders         
weak_cutoff = 2
# maximum percentrank for strong binders
strong_cutoff = 0.5
# optional tsv output name
tsv_output = 'mhc.tsv'
# if to make a tsv summary
make_tsv= True


class Peptideresult:
    '''
    Simply holds all peptide data in one place
    Useful for further annotation or plotting
    '''
    def __init__(self,seq):
        '''
        percentranks - dictionary containing allele:percentrank values
        affinities - dictionary containing allele:affinity values
        sequence -  peptide sequence
        rankcolumn,affinitycolumn - holds values of column numbers holding above values in text file input

        '''
        self.percentranks = {}
        self.affinities = {}
        self.sequence = seq
        self.rankcolumn = 12
        self.affinitycolumn = 15
    def eatline(self, line):
        '''
        adds information from single line

        '''
        line = ' '.join(line.split())
        line = line.split(' ')
        self.percentranks[line[1]] =float(line[self.rankcolumn])
        self.affinities[line[1]] = float(line[self.affinitycolumn])
        self.bestrank = min([self.percentranks[x] for x in self.percentranks.keys()])
    
           

def get_peptidelist(inputdata):
    '''
    generate peptide list from lines of file
    '''
    ifbreak = 0
    peptidelist = []
    for item in inputdata:
        if item[0:3] == "   ":
            line = item
            line = ' '.join(line.split())
            line = line.split(" ")
            pepseq = line[2]
            peptidelist.append(pepseq)
        if item[0] == "-":
                ifbreak  +=1     
        if ifbreak == 3:
            break
    return peptidelist

def eatdata(inputdata, peptidelist):
    '''
    generate list of PeptideResult objects from lines of file and peptide list
    '''
    itemlist = []

    for peptide in peptidelist:
        pep = Peptideresult(peptide)
        inputdata_part = [x for x in inputdata if peptide in x]
        for line in inputdata_part:
            linedata = line.split()
            if len(linedata)>3 and peptide == linedata[2]:
                pep.eatline(line)
        itemlist.append(pep)     
       
    return itemlist



def count_allele(pepdata,allele,minrank_strong,minrank_weak):
    '''
    generate count of weak and stromg binders for allele
    strange string output format is used to make drawing a stacked bar graph easier

    '''
    weakcount = 0
    strongcount = 0
    for peptide in pepdata:
        if peptide.percentranks[allele] <= minrank_weak:
            weakcount+=1
        if peptide.percentranks[allele] <= minrank_strong:
            strongcount+=1
    colval = str(weakcount)+"_"+str(strongcount)
            
        
    return colval
        
def make_netmhcpan_bargraph(dictdata,alleles):
    '''
    stacked bar graph 
    '''
    stack = ["Total Quantified","Binders","Strong Binders"]
    labelvals = ["Total Quantified","Binders","Strong Binders"]
    values_stack = [dictdata[x] for x in stack]
    colors = ["#3C5488","#FA8366","#DC0000"]
    colorvals = ["#3C5488","#FA8366","#DC0000"]
    positions = [1,1,1]
    sizes = [0.5,0.6,0.7]
    a=1
    for allele in alleles:
        a+=1
        positions.append(a)
        colors.append("#FA8366")
        stack.append(allele)
        values_stack.append(int(dictdata[allele].split("_")[0]))
        sizes.append(0.5)
        values_stack.append(int(dictdata[allele].split("_")[1]))
        sizes.append(0.5)
        positions.append(a)
        colors.append("#DC0000")
        stack.append(allele)
        
        
    to_df = {"stack":stack,"values":values_stack,"Binder type":colors,"positions":positions,"size":sizes}
    df = pd.DataFrame(to_df)
    xlabels = ["Total"] + list(alleles)
    xbreaks = [x for x in range(1,len(alleles)+2)]

    
    
    plot = \
    \
    plt.ggplot(df) +\
    plt.aes(x = "positions",y="values", fill = "Binder type") +\
    plt.geom_col(position = "identity")+\
    plt.scale_fill_manual(values =colorvals,labels=labelvals)+\
    plt.labs(x = "", y = "Number of peptides")+\
    plt.scales.scale_y_continuous(expand = (0,0)) +\
    plt.scales.scale_x_continuous(breaks = xbreaks,labels = xlabels ) +\
    plt.theme(panel_background= plt.element_rect(fill = 'white', colour = 'white'),\
    panel_grid= plt.element_line(colour = "#C4C8DA",size =0.5),\
    axis_line_x = plt.element_line(colour = "black" ),\
    axis_line_y = plt.element_line(colour = "black" ),\
    axis_title_y=plt.element_text(size = 14),\
    axis_text_x=plt.element_text(size=10,rotation=60,color="black"),\
    axis_text_y=plt.element_text(size=12,color="black"))

    plot.draw(show =True)
    plot.save('mhc',dpi=1000)
    
def make_dataframe(peptides,outname):
    dataf_source = []
    for peptide in peptides:
        line = {}
        line['Peptide'] = peptide.sequence
        for key in peptide.percentranks.keys():
            line[key+'rank'] = peptide.percentranks[key]
            line[key+'affinity'] = peptide.affinities[key]
        dataf_source.append(line)
    dataf = pd.DataFrame(dataf_source)
    dataf.to_csv(outname,sep='\t',index=False)
        
    
    

allpeps    = []
predictions = []
for item in netmhcpans:
    netmhcpan = open(item, "r").readlines()
    peps = get_peptidelist(netmhcpan)
    preds = eatdata(netmhcpan,peps)
    allpeps = allpeps+peps
    predictions = predictions+preds
strongbinders = []
allbinders = []
for peptide in predictions:
    if peptide.bestrank<=0.5:
        strongbinders.append(peptide)
    if peptide.bestrank<=2:
        allbinders.append(peptide)
        
# prepare data for plotting       
graphdata = {"Total Quantified": total_quantified,"Binders":len(allbinders),"Strong Binders":len(strongbinders)}
alleles = predictions[0].percentranks.keys()
for allele in alleles:
    graphdata[allele] = count_allele(predictions,allele,strong_cutoff,weak_cutoff)
print('Total peptide number ' + str(len(predictions)))    
make_netmhcpan_bargraph(graphdata,alleles)
if make_tsv:
    make_dataframe(predictions,tsv_output)










