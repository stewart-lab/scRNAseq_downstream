# parse DE marker files and create table of cluster, annotation, and marker genes

import os
import sys
import pandas as pd
import operator

def main():
    DIR= str(sys.argv[1]) # directory for DE gene output files
    df0 = pd.DataFrame()
    clust_list = []
    
    for file in os.listdir(DIR): 
        if file.startswith("KnownDE.markers"):
            clust = file.split(".")[1]
            clust = clust.split("_")[-1]
            clust_list.append(clust)
            df = pd.read_csv(DIR+file, sep="\t")
            df = df.sort_values('rank.logFC.cohen')
            #df2 = df.loc[(df['rank.logFC.cohen']<11)]
            celltypelist = df['Cell.type'].unique()
            celltype = "-".join(celltypelist)
            df2 = df.loc[:, ["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen","mean.AUC",
                             "median.AUC","mean.logFC.detected","median.logFC.detected","Cell.type"]]
            for file2 in os.listdir(DIR):
                if file2.startswith("Top100DEgenes_"):
                    clust2 = file2.split(".")[0]
                    clust2 = clust2.split("_")[-1]
                    if clust2 == clust:
                        df3 = pd.read_csv(DIR+file2, sep="\t", header=0) #names=col_names)
                        cols= list(df3.columns)
                        cols.insert(0, "Row.names")
                        df3 = pd.read_csv(DIR+file2, sep="\t", header=0, names=cols)
                        #df3.columns = cols
                        print(df3.columns)
                        df3 = df3.loc[:, ["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen","mean.AUC",
                                          "median.AUC","mean.logFC.detected","median.logFC.detected"]]
                        df3 = df3.loc[(df3['median.logFC.cohen']>=0.35)]
                        df4 = pd.merge(df2, df3, on=["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen",
                                                     "mean.AUC","median.AUC","mean.logFC.detected","median.logFC.detected"], how="outer")
                        print(df4.head())
            
            df3['Cluster']=clust
            df3['Cluster.Cell.type']=celltype
            # df3 = df3.drop(columns=["mean.logFC.cohen_x", "median.logFC.cohen_x", "rank.logFC.cohen_x", "mean.AUC_x", "median.AUC_x", "mean.logFC.detected_x", "median.logFC.detected_x"])
            # df3.rename(columns={'mean.logFC.cohen_y': 'mean.logFC.cohen', 'median.logFC.cohen_y': 'median.logFC.cohen', 
            #                     "rank.logFC.cohen_y": "rank.logFC.cohen", "mean.AUC_y": "mean.AUC", "median.AUC_y": "median.AUC", 
            #                     "mean.logFC.detected_y":"mean.logFC.detected", "median.logFC.detected_y": "median.logFC.detected"}, inplace=True)
            print(df3.head())
            df0 = pd.concat([df0, df3], ignore_index=True)
    for file3 in os.listdir(DIR):
        if file3.startswith("Top100DEgenes_"):
            clust3 = file3.split(".")[0]
            clust3 = clust3.split("_")[-1]
            if clust3 not in clust_list:
                df5 = pd.read_csv(DIR+file3, sep="\t", header=0)
                cols= list(df5.columns)
                cols.insert(0, "Row.names")
                df5 = pd.read_csv(DIR+file3, sep="\t", header=0, names=cols)
                df5 = df5.loc[:, ["Row.names", "mean.logFC.cohen", "median.logFC.cohen", "rank.logFC.cohen","mean.AUC","median.AUC","mean.logFC.detected","median.logFC.detected"]]
                df5 = df5.loc[(df5['median.logFC.cohen']>=0.35)]
                df5['Cluster']=clust3
                df5['Cluster.Cell.type']="unknown"
                df0 = pd.concat([df0, df5], ignore_index=True)
            else:
                pass
    
    df0.to_csv(DIR+"all_clusters_allDEgenes_logFC_0.35.txt", sep="\t", index=False)


if __name__ == '__main__':
    main()