import os, sys
import pandas as pd
startdir = "/w5home/bmoore/Pierre_sc_zebrafish/cellchat2/"
startdir2 = "/w5home/bmoore/Pierre_sc_zebrafish/nichnet_2wk_R-NR_20251110_113947/"
target = ["5","7"]
for dir in os.listdir(startdir):
    if os.path.isdir(os.path.join(startdir, dir)):
        #print(dir)
        for file in os.listdir(startdir + dir):
            if file.endswith("_ligand-recept.txt"):
                print(file)
                filename = startdir + dir + "/" + file
                df = pd.read_csv(filename, sep='\t', header=(0))
                for t in target:
                    df2 = df[df["target"] == t]
                    #df2 = df2.loc[:,['source','target','ligand','receptor']]
                    print(df2)
                    df2['Combined'] = df2['ligand'] + "-" + df2['receptor']
                    LR_result = set(list(df2['Combined']))
                    # lig_result = set(list(df2["ligand"]))
                    # rec_result = set(list(df2["receptor"]))
                    # combine
                    # fin_result = set(lig_result + rec_result)
                
                    # fin_result = pd.DataFrame(data=fin_result, columns=["gene"])
                    # print(fin_result)
                    # output = filename + "_" + t + "_genes.txt"
                    # fin_result.to_csv(output, sep="\t")
                    if t == "5":
                        t2 = "Fibroblasts"
                    elif t == "7":
                        t2 = "T_Cells"
                    else:
                        print(t)
                    for file in os.listdir(startdir2):
                        if file.endswith("_ligand-receptor_links.txt"):
                            print(file)
                            if t2 in file:
                                filename2 = startdir2 + file
                                print(filename2)
                                df3 = pd.read_csv(filename2, sep='\t', header=(0))
                                print(df3)
                                df3['Combined'] = df3['from'] + "-" + df3['to']
                                LR_result2 = set(list(df3["Combined"]))
                                # lig_result2 = set(list(df2["from"]))
                                # rec_result2 = set(list(df2["to"]))
                                # print(lig_result2)
                                # print(rec_result2)
                                fin_LR_result = []
                                for l in LR_result:
                                    if l in LR_result2:
                                        print(l)
                                        fin_LR_result.append(l)
                                fin_LR_result = pd.DataFrame(data=fin_LR_result, columns=["ligand-receptor"])
                                print(fin_LR_result)
                                output = filename + "_" + t + "_LR.txt"
                                fin_LR_result.to_csv(output, sep="\t")
                                