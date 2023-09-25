# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 16:58:23 2022

@author: Jana
"""



# package import
import AccCalc_package_1022 as acc

# genotype in HapMap format upload
file = open("Soy50k_old_hmp.txt","rt")
tab = acc.imp_hmp_tab(file)
file.close()
#print(tab)

# phenotype input
file=open("Flow_col_for_GAPIT.txt","rt")
phen = acc.imp_phen_tab(file,["1"],["2"],Head=True,sep="\t")
#print(phen)

# genotype and phenotype ordering
phen=acc.sort_names(tab,phen)

# Accuracy calculation
Accur = acc.acc_cal(tab, phen, cal="B", flip="on")

# p-value input
file_pval = open("GAPIT.MLM.Flow_col.GWAS.Results.csv", "rt")
Accur = acc.add_pval(Accur,file_pval,file_sep=",",pval_pos=3,chrom_pos=1,pos_pos=2)

# Accu tab rounding (only for space saving an more easy to read results)
Accur=Accur.round({"p-value":8,"-log10pval":2,"Avr_acc(%)":2,"Acc_WT(%)":2,"Acc_MUT(%)":2,
          "Avr_acc_pes(%)":2,"Acc_pes_WT(%)":2,"Acc_pes_MUT(%)":2,"Comb_acc(%)":2,
          "Comb_acc_pes(%)":2,"Missing_gen(%)":2,"Missing_gen_for_WT(%)":2,
         " Missing_gen_for_MUT_(%)":2,"Missing_phen(%)":2}) 

# adding distance in bp to the selected position
Accur=acc.add_distance(Accur,point="max_pval")
Accur.to_csv("Acc_Soy50kOld_Flow_col.csv") # saving the full results

# visualization
acc.plot_accuracy(Accur)  

# selecting and ordering the results for only the important ones (hardly dependent on data)
Accur.to_csv("Acc_Soy50kOld_Flow_col_markers_all.csv", index=False)
Accur = Accur[Accur["Avr_acc(%)"]>80]
Accur = Accur.sort_values(by=["-log10pval",'Avr_acc(%)'],ascending=False)

# saving the selected results
Accur.to_csv("Acc_Soy50kOld_Flow_col_markers_ACC80.csv", index=False)



