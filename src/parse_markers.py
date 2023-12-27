# parse DE marker files and create table of cluster, annotation, and marker genes

import os
import sys
import pandas as pd
import operator

def main():
    DIR= "/w5home/bmoore/scRNAseq/GAMM/GAMM_S2/output_20230830_155530/"#sys.argv[1]
    df0 = pd.DataFrame()
    clust_list = []
    for file in os.listdir(DIR): 
        if file.startswith("gamm.known"):
            clust = file.split(".")[2]
            clust = clust.split("_")[-1]
            clust_list.append(clust)
            df = pd.read_csv(DIR+file, sep="\t")
            df = df.sort_values('rank.logFC.cohen')
            df2 = df.loc[(df['rank.logFC.cohen']<11)]
            celltypelist = df2['Cell.type'].unique()
            celltype = "-".join(celltypelist)
            df3 = df.loc[:, ["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen","mean.AUC",
                             "median.AUC","mean.logFC.detected","median.logFC.detected","Cell.type"]]
            for file2 in os.listdir(DIR):
                if file2.startswith("gamm.top100genes"):
                    clust2 = file2.split(".")[1]
                    clust2 = clust2.split("_")[-1]
                    if clust2 == clust:
                        col_names=["Row.names", "self.average", "other.average", "self.detected", "other.detected", 
                                   "mean.logFC.cohen", "min.logFC.cohen", "median.logFC.cohen", "max.logFC.cohen", 
                                   "rank.logFC.cohen", "full.logFC.cohen.X2", "full.logFC.cohen.X3", "full.logFC.cohen.X4", 
                                   "full.logFC.cohen.X5", "full.logFC.cohen.X6", "full.logFC.cohen.X7", "full.logFC.cohen.X8", 
                                   "full.logFC.cohen.X9", "full.logFC.cohen.X10", "full.logFC.cohen.X11", "full.logFC.cohen.X12", 
                                   "full.logFC.cohen.X13", "full.logFC.cohen.X14", "full.logFC.cohen.X15", "mean.AUC", 
                                   "min.AUC", "median.AUC", "max.AUC", "rank.AUC", "full.AUC.X2", "full.AUC.X3", "full.AUC.X4", 
                                   "full.AUC.X5", "full.AUC.X6", "full.AUC.X7", "full.AUC.X8", "full.AUC.X9", "full.AUC.X10", 
                                   "full.AUC.X11", "full.AUC.X12", "full.AUC.X13", "full.AUC.X14", "full.AUC.X15", 
                                   "mean.logFC.detected", "min.logFC.detected", "median.logFC.detected", "max.logFC.detected", 
                                   "rank.logFC.detected", "full.logFC.detected.X2", "full.logFC.detected.X3", "full.logFC.detected.X4", 
                                   "full.logFC.detected.X5", "full.logFC.detected.X6", "full.logFC.detected.X7", "full.logFC.detected.X8", 
                                   "full.logFC.detected.X9", "full.logFC.detected.X10", "full.logFC.detected.X11", 
                                   "full.logFC.detected.X12", "full.logFC.detected.X13", "full.logFC.detected.X14", "full.logFC.detected.X15"]
                        df4 = pd.read_csv(DIR+file2, sep="\t", header=0, names=col_names)
                        df4 = df4.loc[:, ["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen","mean.AUC","median.AUC","mean.logFC.detected","median.logFC.detected"]]
                        df4 = df4.loc[(df4['median.logFC.cohen']>=0.35)]
                        df3 = pd.merge(df3, df4, on=["Row.names"], how="outer")
                        print(df3.head())
            
            df3['Cluster']=clust
            df3['Cluster.Cell.type']=celltype
            df3 = df3.drop(columns=["mean.logFC.cohen_x", "median.logFC.cohen_x", "rank.logFC.cohen_x", "mean.AUC_x", "median.AUC_x", "mean.logFC.detected_x", "median.logFC.detected_x"])
            df3.rename(columns={'mean.logFC.cohen_y': 'mean.logFC.cohen', 'median.logFC.cohen_y': 'median.logFC.cohen', 
                                "rank.logFC.cohen_y": "rank.logFC.cohen", "mean.AUC_y": "mean.AUC", "median.AUC_y": "median.AUC", 
                                "mean.logFC.detected_y":"mean.logFC.detected", "median.logFC.detected_y": "median.logFC.detected"}, inplace=True)
            print(df3.head())
            df0 = pd.concat([df0, df3], ignore_index=True)
    for file3 in os.listdir(DIR):
        if file3.startswith("gamm.top100genes"):
            clust3 = file3.split(".")[1]
            clust3 = clust3.split("_")[-1]
            if clust3 not in clust_list:
                df5 = pd.read_csv(DIR+file3, sep="\t", header=0, names=col_names)
                df5 = df5.loc[:, ["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen","mean.AUC","median.AUC","mean.logFC.detected","median.logFC.detected"]]
                df5 = df5.loc[(df5['median.logFC.cohen']>=0.5)]
                df5['Cluster']=clust3
                df5['Cluster.Cell.type']="unknown"
                df0 = pd.concat([df0, df5], ignore_index=True)
            else:
                pass
    df0.to_csv(DIR+"gamm.all_clusters_allgenes.txt", sep="\t", index=False)

if __name__ == '__main__':
    main()