import os, sys
import pandas as pd
startdir = "/w5home/bmoore/Pierre_sc_zebrafish/cellchat2/"
startdir2 = "/w5home/bmoore/Pierre_sc_zebrafish/nichnet_4wk_R-NR_20251215_151657/"
timepoint = "4Wk"
target = ["5","7"]
gpt_file= "extracted_kmGPT_rawdata_allweeks.csv"

# get gpt df
gpt_df = pd.read_csv(startdir + gpt_file, sep=',', header=(0))
# set final rsult to 0
result2 = pd.DataFrame()

# loop trhough files
for file in os.listdir(startdir):
    if os.path.isfile(os.path.join(startdir, file)):
        if "cellchat_df_net_ligand-recept_" and timepoint in file:
            print(file)
            if "WkNR" in file:
                condition = "NR"
            elif "WkR" in file:
                condition = "R"
            else:
                print(file, "no condition")
            if "Fibroblast" in file:
                t = "7"
                t2 = "Fibroblasts"
                t3 = "Fibroblast"
            elif "Tcell" in  file:
                t = "5"
                t2 = "T_Cells"
                t3 = "Tcell"
            else:
                print("no target info")
            filename = startdir + file
            
            # get cell chat result
            df = pd.read_csv(filename, sep='\t', header=(0))
            print(df)
            # subset by target (should already be subsetted)
            df["target"] = df["target"].astype(str)
            df2 = df[df["target"] == t]
            print("subsetted ",df2)
            # get ligand and receptr from cell chat
            df2['Combined'] = df2['ligand'] + "-" + df2['receptor']
            LR_result = set(list(df2['Combined']))
                
            # get niche net results
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
                        LR_result3 = []
                        for l in LR_result:
                            if l in LR_result2:
                                print(l)
                                LR_result3.append(l)
                                
                        # subset gpt file
                        print("gpt df: ", gpt_df)
                        gpt_df_sub = gpt_df[(gpt_df["celltype"] == t3) & (gpt_df["timepoint"] == timepoint.lower()) & (gpt_df["score"] > 0)]
                        print("gpt sub df: ",gpt_df_sub, set(gpt_df_sub["b_term"]))
                        fin_LR_result = {}
                        for lr in LR_result3:
                            lig = lr.split("-")[0].strip()
                            rec = lr.split("-")[1].strip()
                            
                            
                            if str(lig) in set(gpt_df_sub["b_term"]):
                                print(lig)
                                print(lr)
                                index_list = gpt_df_sub.isin([str(lig)]).any(axis=1).index.tolist()
                                print(index_list)
                                a_term_list = []
                                for index in index_list:
                                    a_term = gpt_df_sub["a_term"][index]
                                    if a_term not in a_term_list:
                                        a_term_list.append(a_term)
                                print(a_term_list)
                                a_term_str = ",".join(a_term_list)
                                fin_LR_result[lr] = a_term_str
                            elif str(rec) in set(gpt_df_sub["b_term"]):
                                print(rec)
                                print(lr)
                                index_list = gpt_df_sub.isin([str(lig)]).any(axis=1).index.tolist()
                                print(index_list)
                                a_term_list = []
                                for index in index_list:
                                    a_term = gpt_df_sub["a_term"][index]
                                    if a_term not in a_term_list:
                                        a_term_list.append(a_term)
                                print(a_term_list)
                                a_term_str = ",".join(a_term_list)
                                fin_LR_result[lr] = a_term_str
                            else:
                                print("lig-rec not in gpt df, ", lr, str(lig), str(rec))
                                #pass
                        fin_LR_result = pd.DataFrame.from_dict(data=fin_LR_result, orient='index', columns=["a_terms"])
                        fin_LR_result["Condition"]= condition
                        fin_LR_result["Timepoint"]= timepoint 
                        fin_LR_result["Target"]= t3
                        print(fin_LR_result)
                        # change index to col name
                        #fin_LR_result = fin_LR_result.set_index('ligand-receptor')
                        #fin_LR_result['ligand-receptor'] = fin_LR_result.index
                        # reset
                        fin_LR_result = fin_LR_result.reset_index(level=0, names='ligand-receptor')
                        print(fin_LR_result)
                        print(df2)
                        result = pd.merge(fin_LR_result, df2, how="left", left_on='ligand-receptor', right_on="Combined", sort=True) #left_index=True
                        print(result)
                        if not result2.empty:
                            result2 = pd.concat([result2, result], ignore_index=True)
                        else:
                            result2 =result
                            
                            
                            
output = "/w5home/bmoore/Pierre_sc_zebrafish/cellchat-nichenet-kmgpt_" + timepoint + "_LR.txt"
result2.to_csv(output, sep="\t")
                                