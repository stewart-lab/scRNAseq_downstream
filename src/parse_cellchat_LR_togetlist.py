import os, sys
import pandas as pd
startdir = "/w5home/bmoore/Pierre_sc_zebrafish/cellchat2/"
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
                    lig_result = list(df2["ligand"])
                    rec_result = list(df2["receptor"])
                    # combine
                    fin_result = set(lig_result + rec_result)
                
                    fin_result = pd.DataFrame(data=fin_result, columns=["gene"])
                    print(fin_result)
                    output = filename + "_" + t + "_genes.txt"
                    fin_result.to_csv(output, sep="\t")