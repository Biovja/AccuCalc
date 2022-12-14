# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 10:50:22 2022

@author: Jana
"""

## An abbreviations used in the code notes
#### Ref - reference allele
#### Alt - alternate allele
#### WT - wild type phenotype
#### MUT - mutant phenotype
#### Acc - accuracy


#Import packages 
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import StrMethodFormatter

################################################################################################
################################################################################################
def vcfline_to_tab (line,alt,het):
    try: int(line[0])
    except:
        ch=""
        for v in line[0]:
            if v.isdigit(): ch=ch+v
        line[0]=ch
    if alt!="max":
            ln = line[:9]
            for a in line[9:]:
                if a[0]=="0" and a[2]=="0": ln.append(0)
                elif a[0]=="1" and a[2]=="1": ln.append(1)
                elif alt=="all" and het=="N":
                    try: 
                        if int(a[0])==int(a[2]): ln.append(1)
                        else: ln.append(-2)
                    except: ln.append(-2)
                elif alt=="all":
                    try: 
                        if int(a[0])==int(a[2]): ln.append(1)
                        elif het=="Ref": ln.append(0)
                        elif het=="Alt": ln.append(1)
                    except: ln.append(-2)
                else: 
                    if het!="N":
                        try: 
                         int(a[0]),int(a[2])
                         if het=="Ref": ln.append(0)
                         elif het=="Alt": ln.append(1)
                        except: ln.append(-2)
                    else: ln.append(-2)
    elif alt=="max":
            ln = line[:9]
            d=dict()
            for a in line[9:]:
                try: 
                    if int(a[0])==int(a[2]) and a[0]!="0": d[a[0]] = d.get(a[0],0)+1
                except: continue
            if len(d)>=1: sel_alt = max(d, key=d.get)
            else: sel_alt=1
            for a in line[9:]:
                if a[0]=="0" and a[2]=="0": ln.append(0)
                elif a[0]==sel_alt and a[2]==sel_alt: ln.append(1)
                else: 
                    if het!="N":
                         if het=="Ref" and (a[0]==sel_alt or a[2]==sel_alt): ln.append(0)
                         elif het=="Alt" and (a[0]==sel_alt or a[2]==sel_alt): ln.append(1)
                         else: ln.append(-2)
                    else: ln.append(-2)
    return ln

##################################################################################################
def hmpline_to_tab (line):
    allele=line[1].split("/")
    try: allele[1]
    except: allele.append("NA")
    inf="strand="+line[4]+";assembly="+line[5]+";center="+line[6]+";protLSID="+line[7]+";assayLSID="+line[8]+";panelLSID="+line[9]+";QCcode="+line[10]
    ln=[line[2],line[3],line[0],allele[0],allele[1],"NA","NA",inf,"NA"]
    for a in line[11:]:
        if a==allele[0] or a==allele[0]+allele[0]: ln.append(0)
        elif a!="N" and a!="NA" and a!="NaN" and a!="*" and a!="." and a!="-": ln.append(1)
        else: ln.append(-2)
    return ln
##################################################################################################
## Import hapmap file into a numerically coded table in Pandas (DataFrame),  
### table format suitable for the other functions in the package
######## 0 for the first allele in hmp, 1 for the second allele, -2 for missing and heterozygotes
## area: allow to select only a small part from the vcf file
#### list of integers with a genomic area definition in format: [chromosome,start position,stop position]
def imp_hmp_tab(file, area="off"): 
    n,lst=0,list()
    for  line in file:
        line=line.split()
        
        if n==1:
            if area!="off" and int(line[2])==area[0] and int(line[3])>area[1] and int(line[3])<area[2]:
                lst.append(hmpline_to_tab (line))    
            elif area=="off":
                lst.append(hmpline_to_tab (line))
                
        else:
            hd=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+line[11:]
            n=1
    
    gen_tab=pd.DataFrame(lst,columns=hd)   
    gen_tab["#CHROM"]=pd.to_numeric(gen_tab["#CHROM"]) 
    gen_tab["POS"]=pd.to_numeric(gen_tab["POS"]) 
    return gen_tab
    
##################################################################################################
## Import vcf file into a numericaly coded table in Pandas (DataFrame),  
### table format suitable for the other functions in the package
######## 0 for Ref (0/0), 1 for Alt, -2 for missing and heterozygotes
## area: allow to select only a small part from the vcf file
#### list of integers with a genomic area definition in format: [chromosome,start position,stop position]
## alt: define Alt coding in case of more than one Alt variant
#### all - all Alt's coded as Alt(1), max - only the most frequent Alt as Alt, 1st - only the first (1/1) Alt as Alt
## het: define coding for heterozygosity
#### N - defined as missing, Ref - defined as reference, Alt - defined as Alt
def imp_vcf_tab(file, area="off", alt="all", het="N"):
    lst = list()
    start = "n"
    for line in file:
        if line.startswith("#CHROM"):
            Head = line.split()
            start="y"
        elif start=="y":
            
            line=line.split()
            if area=="off": lst.append(vcfline_to_tab (line,alt,het))
            else:
                try:
                    if int(line[0])==area[0] and int(line[1])>area[1] and int(line[1])<area[2]:
                        lst.append(vcfline_to_tab (line,alt))
                        
                except: 
                    ch=""
                    
                    for v in line[0]:
                        if v.isdigit(): ch=ch+v
                    if ch!="":
                        
                        if int(ch)==area[0] and int(line[1])>area[1] and int(line[1])<area[2]:
                                line[0]=ch
                                lst.append(vcfline_to_tab (line,alt,het))
               
    gen_tab=pd.DataFrame(lst,columns=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+Head[9:])   
    gen_tab["#CHROM"]=pd.to_numeric(gen_tab["#CHROM"]) 
    gen_tab["POS"]=pd.to_numeric(gen_tab["POS"]) 
    return gen_tab
##################################################################################################
## Import phenotype file into a numerically coded table in Pandas (DataFrame),  
### table format suitable for the other functions in the package
######## 0 for WT, 1 for MUT, -2 for missing and other phenotypes
## wt/mut: list of phenotypes used as WT/MUT 
### ect. ["B","Br"] for categorical or [1,10] for quantitative (numbers sets a range)
## phen_type: determinates phenotype type 
### "categ" for categorical or binary phenotype or "quant" for quantitative phenotype
## Head=True/False: default True 
## sep: type of separator in phenotype file
## col_n: column number in phenotype file, counted from 0 (sample names in the first colum)
def imp_phen_tab(file,wt,mut,phen_type="categ",Head=True,sep="\t",col_n=1):
    phn,hd,n=list(),list(),0
    for line in file:
        line=line.split(sep)
        if Head and n==0:n=1
        elif phen_type=="quant":
            hd.append(line[0])
            try:
                if float(line[col_n])>=wt[0] and float(line[col_n])<wt[1]: phn.append(0)
                elif float(line[col_n])>=mut[0] and float(line[col_n])<mut[1]: phn.append(1)
                else: phn.append(-2)
            except: phn.append(-2)
        else:
            hd.append(line[0])
            
            if line[col_n].strip() in wt: phn.append(0)
            elif line[col_n].strip() in mut: phn.append(1)
            else: phn.append(-2)
    phen_tab=pd.DataFrame([phn],columns=hd)
    return phen_tab
##################################################################################################
## Create Synthetic phenotype from selected genotype line in vcf file
## return Phenotype table with numerically coded phenotype info (DataFrame),  
### table format suitable for the other functions in the package
######## 0 for WT, 1 for MUT, -2 for missing and other phenotypes
## chrom/pos numeric specification for line selection
## wt/mut: list of phenotypes used as WT/MUT etc. ["0/0","0/1"]
def synthetic_phen_from_vcf(file,chrom,pos,wt,mut):
    start,phn,hd="n",list(),list()
    for line in file:
        if start=="y":
            line=line.split()
            if line[0].isdigit():
                if chrom==int(line[0]) and pos==int(line[1]):
                    for a in line[9:]:
                        if a[:3] in wt:phn.append(0)
                        elif a[:3] in mut: phn.append(1)
                        else: phn.append(-2)
                    break
            else:
                ch=""
                for v in line[0]:
                    if v.isdigit(): ch=ch+v
                if chrom!="" and chrom==int(ch) and pos==int(line[1]):
                    for a in line[9:]:
                        print(a[:3])
                        if a[:3] in wt:phn.append(0)
                        elif a[:3] in mut: phn.append(1)
                        else: phn.append(-2)
                    break         
        elif line.startswith("#C"): start,hd="y",line.split()[9:]
    
    if phn==[]: return print("ERROR! Data for synthetic phenotype was not found.")
    else:
        phen_tab = pd.DataFrame([phn],columns=hd)
        return phen_tab
##################################################################################################
## Import phenotype file into a numerically coded table in Pandas (DataFrame),  
### table format suitable for the other functions in the package
######## 0 for WT, 1 for MUT, -2 for missing and other phenotypes
## chrom: integer value of chromosome of selected locus
## pos: integer value of the position of selected locus
## wt/mut: integer value (0 for Ref and 1 for Alt) of WT/MUT 
def synthetic_phen_from_tab(gen_tab,chrom,pos,wt,mut):
    if wt==0 and mut==1:
        phen_tab = gen_tab[(gen_tab["#CHROM"]==chrom) & (gen_tab["POS"]==pos)]
        phen_tab = (phen_tab.iloc[:, 9:])
    else:
        hd = list(gen_tab.columns.values)
        ln = gen_tab[(gen_tab["#CHROM"]==chrom) & (gen_tab["POS"]==pos)].values.flatten().tolist()
        phn=list()
        for a in ln[9:]:
            if a==wt: phn.append(0)
            elif a==mut: phn.append(1)
            else: phn.append(-2)
        phen_tab = pd.DataFrame([phn],columns=hd[9:])
    return phen_tab
#####################################################################################################
## Write a phenotype txt file from pheno_tab
## name: name of the txt phenotype file
## pheno: Name of phenotype column in txt phenotype file
def export_phen_tab(phen_tab,name,pheno="Phenotype"):
    names=list(phen_tab.columns.values)
    phn = phen_tab.iloc[0].values.flatten().tolist()
    res=open(f"{name}.txt","w")
    res.write(f"Sample\t{pheno}\n")
    for name,ph in zip(names,phn):
        res.write(f"{name}\t{ph}\n")
    return print("Done.")
#####################################################################################################
##Create Acc_tab 
### cal: A - Ref is WT and Alt is MUT, B - Ref is MUT and Alt is WT
### flip: on - if Acc_real < 50% flip model A/B
def acc_cal(gen_tab, phen_tab, cal="A", flip="off"):
    if list(gen_tab.columns.values)[9:]==list(phen_tab.columns.values):
        phen = np.array(phen_tab)
        phenX = 1-phen
        phen_miss=np.count_nonzero(phen==-2)/phen.size*100
        n_wt,n_mut= np.count_nonzero(phen==0),np.count_nonzero(phen==1)
        acc_real,acc_pes,acc_realc,acc_pesc,ref_phn,gen_miss=list(),list(),list(),list(),list(),list()
        n_ref,n_alt,acc_wt_real,acc_mut_real,acc_wt_pes,acc_mut_pes,=list(),list(),list(),list(),list(),list()
        wt_miss,mut_miss=list(),list()
        for n in range(0,len(gen_tab.index.values)):
             gen = np.array(gen_tab.iloc[n][9:])
             l_all = len(gen)
             n_ref.append(np.count_nonzero(gen==0)),n_alt.append(np.count_nonzero(gen==1))
             gen_miss.append((np.count_nonzero(gen==-2)/l_all)*100)
             
             res = gen+phen
             A,B,M = np.count_nonzero(res==0),np.count_nonzero(res==2),np.count_nonzero(res<0)
             resX = gen+phenX
             C,D = np.count_nonzero(resX==0),np.count_nonzero(resX==2)
             ABCD = A+B+C+D
             AD,BC = A+D,B+C
                         
             if l_all!=(A+B+C+D+M): print("ERROR in counts!",l_all,(A+B+C+D+M))
             
             if cal=="A":
                 ref_phn_val="WT"
                 if (AD)!=0 and (BC)!=0 :
                     accrw,accrm=(A/AD)*100,(B/BC)*100
                     accr=(accrw+accrm)/2
                 else: accr,accrw,accrm=0,0,0
                 accpw,accpm=(A/n_wt)*100,(B/n_mut)*100
                 accp=(accpw+accpm)/2
                 wtm,mutm=((n_wt-AD)/n_wt)*100,((n_mut-BC)/n_mut)*100
                 if (ABCD)!=0: accrc=(A+B)/(ABCD)*100
                 else: accrc=0
                 accpc=(A+B)/l_all*100
                 if flip=="on" and accr<50:
                     if (AD)!=0 and (BC)!=0:
                         accrw,accrm=(D/AD)*100,(C/BC)*100
                         accr=(accrw+accrm)/2
                     else: accr,accrw,accrm=0,0,0
                     accpw,accpm=(D/n_wt)*100,(C/n_mut)*100
                     accp=(accpw+accpm)/2
                     if (ABCD)!=0:accrc=(C+D)/(ABCD)*100
                     else: accrc=0
                     accpc=(C+D)/l_all*100
                     ref_phn_val="MUT"
             elif cal=="B":
                 ref_phn_val="MUT"
                 if (AD)!=0 and (BC)!=0:
                     accrw,accrm=(D/AD)*100,(C/BC)*100
                     accr=(accrw+accrm)/2
                 else: accr,accrw,accrm=0,0,0
                 accpw,accpm=(D/n_wt)*100,(C/n_mut)*100
                 accp=(accpw+accpm)/2
                 wtm,mutm=((n_wt-AD)/n_wt)*100,((n_mut-BC)/n_mut)*100
                 if (ABCD)!=0:accrc=(D+C)/(ABCD)*100
                 else: accrc=0
                 accpc=(D+C)/l_all*100
                 if flip=="on" and accr<50:
                     if (AD)!=0 and (BC)!=0:
                         accrw,accrm=(A/AD)*100,(B/BC)*100
                         accr=(accrw+accrm)/2
                     else: accr,accrw,accrm=0,0,0
                     accpw,accpm=(A/n_wt)*100,(B/n_mut)*100
                     accp=(accpw+accpm)/2
                     if (ABCD)!=0:accrc=(A+B)/(ABCD)*100
                     else: accrc=0
                     accpc=(A+B)/l_all*100
                     ref_phn_val="WT"
                          
             acc_real.append(accr),acc_pes.append(accp),acc_pesc.append(accpc),acc_realc.append(accrc),
             ref_phn.append(ref_phn_val),acc_wt_real.append(accrw),acc_mut_real.append(accrm),
             acc_wt_pes.append(accpw),acc_mut_pes.append(accpm),wt_miss.append(wtm),mut_miss.append(mutm)
  
        acc_tab = pd.DataFrame({"#CHROM":gen_tab["#CHROM"],"POS":gen_tab["POS"],"ID":gen_tab["ID"],
                                "REF":gen_tab["REF"],"ALT":gen_tab["ALT"],"INFO":gen_tab["INFO"],
                                'Avr_acc(%)': acc_real,'Acc_WT(%)': acc_wt_real,'Acc_MUT(%)': acc_mut_real,
                                'Avr_acc_pes(%)': acc_pes,'Acc_pes_WT(%)': acc_wt_pes,'Acc_pes_MUT(%)': acc_mut_pes,
                                'Comb_acc(%)': acc_realc, 'Comb_acc_pes(%)': acc_pesc,
                                "Missing_gen(%)" :gen_miss,"Count_ref":n_ref,"Count_alt":n_alt,"Ref_phen_for_Acc":ref_phn,
                                "Missing_gen_for_WT(%)" :wt_miss,"Missing_gen_for_MUT_(%)" :mut_miss,
                                "Missing_phen(%)" :phen_miss,"Count_WT":n_wt,"Count_MUT":n_mut})
        
        acc_tab["#CHROM"]=pd.to_numeric(acc_tab["#CHROM"] ,errors='coerce')
        acc_tab["POS"]=pd.to_numeric(acc_tab["POS"] ,errors='coerce')
        return acc_tab
    else:
        print("Error. Samples in phen_table and gen_table do not fit. Try to use addjust phenotype with sort_names function.")
        print(gen_tab.columns.values[9:])
        print(phen_tab.columns.values)
        print(len(gen_tab.columns.values[9:]),len(phen_tab.columns.values))
#####################################################################################################
## Sort phen_tab by sample order of the gen_tab
### Necessary for Accurycy calculation
def sort_names(gen_tab,phen_tab):
    gen_names=list(gen_tab.columns.values)
    try: 
        phen_tab = phen_tab[gen_names[9:]]
        return phen_tab
    except:
        phen_names=list(phen_tab.columns.values)
        for name in phen_names:
            if name not in gen_names[9:]:
                return print("ERROR. Incorrect sample/sample name in phenotype file. Sample sorting failed!.")
        for name in gen_names:
            if name not in phen_names:phen_tab[name]=[-2]
        phen_tab = phen_tab[gen_names[2:]]
        return phen_tab
    
#####################################################################################################
## Add p-values into Accu_tab
### pval_file: name of file with p-values in format .txt or .csv
### file_sep: type of column separator in pval_file
### pval_pos: position of column with the p-values in pval_file (start from 0)
### chrom_pos: position of column with the chromosome in pval_file (start from 0)
### pos_pos: position of column with the position in pval_file (start from 0)
def add_pval(acc_tab,pval_file,file_sep=",",pval_pos=3,chrom_pos=0,pos_pos=1):
    d_pvl = {}
    for line in pval_file:
        line=line.split(file_sep)
        try: int(line[chrom_pos]),int(line[pos_pos])
        except: continue
        if int(line[chrom_pos]) in d_pvl: d_pvl[int(line[chrom_pos])][int(line[pos_pos])]=line[pval_pos]
        else: d_pvl[int(line[chrom_pos])]={int(line[pos_pos]):line[pval_pos]}
    pval_file.close()
    
    chr_lst=acc_tab["#CHROM"].to_list()
    pos_lst=acc_tab["POS"].to_list()
    pval_list=list()
    for chrom,pos in zip(chr_lst,pos_lst):
        try: pval=d_pvl[chrom][pos]
        except: pval="NA"
        pval_list.append(pval)
       
    acc_tab.insert(6,"p-value",pval_list) 
    acc_tab["p-value"]=pd.to_numeric(pval_list ,errors='coerce')
    acc_tab.insert(7, "-log10pval", pd.to_numeric(-np.log10(acc_tab["p-value"]) ,errors='coerce'))   
    
    return acc_tab

#####################################################################################################   
## Add distance from the position into Accu_tab
### point: position for counting distance
####### max_pval - automatically set as a position with the biggest -log10 pval
####### [chrom,pos] - position defined by the user
def add_distance(acc_tab,point="max_pval"):
    if point=="max_pval":
        try:
            acc_tab = acc_tab.sort_values(by=["-log10pval"],ascending=False)
        except: raise ("ERROR. No p-values available.")
        chrom = acc_tab.iloc[0]["#CHROM"]
        pos = acc_tab.iloc[0]["POS"]
    else:
        try: 
            chrom = int(point[0])
            pos = int(point[1])
        except: raise ("ERROR. Invalid values or format for position assignment.")
    chr_lst=acc_tab["#CHROM"].to_list()
    pos_lst=acc_tab["POS"].to_list()
    
    chr_lst=acc_tab["#CHROM"].to_list()
    pos_lst=acc_tab["POS"].to_list()
    dist_list=list()
    for ch,p in zip(chr_lst,pos_lst):
        if ch!=chrom: dist_list.append("NA")
        else: dist_list.append(abs(p-pos))
    
       
    acc_tab["dist_to_position"]=pd.to_numeric(dist_list ,errors='coerce')
    acc_tab["point"]=f"{chrom},{pos}"
    return acc_tab
#####################################################################################################
# Create predefined plots from accuracy results, before usage add p-values into Acc_tab (function add_pval())
def plot_accuracy(acc_tab):
    chromosomes=np.unique(acc_tab["#CHROM"])
    norm = Normalize(vmin=0, vmax=100)
    if len(chromosomes)>1:
        #print("More than 1 chromosome.")
        acc_tab=acc_tab.sort_values(by=["#CHROM","POS"])
        acc_tab['ind'] = range(len(acc_tab))
        for accuracy in ['Avr_acc(%)','Avr_acc_pes(%)','Comb_acc(%)','Comb_acc_pes(%)']:
            fig = plt.figure(figsize=(9, 3), dpi=300)
            ax = fig.add_subplot()
            cmap = LinearSegmentedColormap.from_list("name",colors=["gainsboro","gainsboro","silver","silver","lightcoral","red","brown"])
            fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),orientation='vertical', label=accuracy,ax=ax)
            cmap = LinearSegmentedColormap.from_list("name",colors=["gainsboro","gainsboro","silver","silver","lightblue","blue","navy"])
            fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),orientation='vertical', label=accuracy,ax=ax)
            
            x_labels,x_labels_pos = list(),list()
            for num,chrom in enumerate(chromosomes):
                df1=acc_tab[acc_tab["#CHROM"]==chrom]
                if num%2==0: cmap = LinearSegmentedColormap.from_list("name",colors=["gainsboro","gainsboro","silver","silver","lightcoral","red","brown"])
                else: cmap = LinearSegmentedColormap.from_list("name",colors=["gainsboro","gainsboro","silver","silver","lightblue","blue","navy"])
                plt.scatter(x=df1['ind'], y=df1['-log10pval'],c=df1[accuracy], cmap=cmap,norm=norm,)
                #group.plot(kind='scatter', x='ind', y='pval',z="Acc",color=cmap, marker=marker[num % len(marker)], ax=ax)
                x_labels.append(chrom)
                x_labels_pos.append((df1['ind'].iloc[-1] + df1['ind'].iloc[0])/2)
            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels)
            plt.title("Accuracy Plot - Genomewise")
            plt.ylabel("-log10(p-val)")
            ax.set_xlabel('Chromosome no.')
            plt.savefig(f"{accuracy}_plot_genomewise.png", dpi=300, bbox_inches='tight')
            plt.show()        
            
            for num,chrom in enumerate(chromosomes):
                df1=acc_tab[acc_tab["#CHROM"]==chrom]
                fig, ax = plt.subplots()
                fig.autofmt_xdate()
                cmap = LinearSegmentedColormap.from_list("name",colors=["lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","grey","grey","grey","grey","grey","darkviolet","blue","c","limegreen","gold","r"])
                fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),orientation='vertical', label=accuracy,ax=ax)
                plt.scatter(df1['POS'], df1['-log10pval'], c=df1[accuracy], alpha=0.7, cmap=cmap,norm=norm)
                ax.set_xticklabels(ax.get_xticks(), rotation = 45)
                ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) 
                plt.title(f"Accuracy Plot - Chromosome {chrom}")
                plt.xlabel(f"Position on chromosome {chrom} (bp)")
                plt.ylabel("-log10(p-val)")
                plt.savefig(f"{accuracy}_plot_chromosome{chrom}.png", dpi=300, bbox_inches='tight')
                plt.show() 
    else:
        for accuracy in ['Avr_acc(%)','Avr_acc_pes(%)','Comb_acc(%)','Comb_acc_pes(%)']:
            try:chrom = chromosomes[0]
            except: chrom = acc_tab["#CHROM"][0]
            fig, ax = plt.subplots()
            fig.autofmt_xdate()
            cmap = LinearSegmentedColormap.from_list("name",colors=["lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","grey","grey","grey","grey","grey","darkviolet","blue","c","limegreen","gold","r"])
            fig.colorbar(ScalarMappable(norm=norm, cmap=cmap),orientation='vertical', label=accuracy,ax=ax)
            plt.scatter(acc_tab['POS'], acc_tab['-log10pval'], c=acc_tab[accuracy], alpha=0.7, cmap=cmap,norm=norm)
            ax.set_xticklabels(ax.get_xticks(), rotation = 45)
            ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) 
            plt.title(f"Accuracy Plot - Chromosome {chrom}")
            plt.xlabel(f"Position on chromosome {chrom} (bp)")
            plt.ylabel("-log10(p-val)")
            plt.savefig(f"{accuracy}_plot_chromosome{chrom}.png", dpi=300, bbox_inches='tight')
            plt.show() 
            