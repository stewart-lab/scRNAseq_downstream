#This script joins 2 data frames together. -df1= first dataframe; -df2 second dataframe; 
#-type {i= inner/intersection, o= outer/union, df1= exact join from first dataframe
# for type inner (i), assumes 2 data frames with class in column 1 for each. for o OR df1, assumes class column1 for df1, but not for df2
#-lower= T or true if you want to adjust dataframes so indexes are all lowercase
import sys, os
import pandas as pd

DF3 = "none"
LOWER= "F"
for i in range (1,len(sys.argv),2):
            if sys.argv[i] == "-df1":
                DF1 = sys.argv[i+1]
            if sys.argv[i] == "-df2":
                DF2 = sys.argv[i+1]
            if sys.argv[i] == "-df3":
                DF3 = sys.argv[i+1]
            if sys.argv[i] == "-type":
                TYPE = str(sys.argv[i+1])
            if sys.argv[i] == "-lower":
                LOWER = str(sys.argv[i+1])
df1 = pd.read_csv(DF1, sep=' ', index_col = 0)
print(list(df1.index))
#df1.drop_duplicates(inplace=True)
df2 = pd.read_csv(DF2, sep='\t', index_col = 0)
print(list(df2.index))
# drop nas
df2.dropna(how='any', axis=0, inplace=True)
# drop duplicates
#df2.drop_duplicates(subset='Mouse.gene.name',inplace=True)
# set index
#df2.set_index('Mouse.gene.name',inplace=True)
#print(list(df2.index))

df1_name = os.path.basename(DF1)
df2_name = os.path.basename(DF2)
#df2.drop_duplicates(inplace=True)
col_list1= list(df1.columns.values)
col_list2= list(df2.columns.values)
x= len(col_list1)
y= len(col_list2)
if LOWER == "T":
	df1.index = df1.index.str.lower()
	df2.index = df2.index.str.lower()
elif LOWER.lower() == "true":
	df1.index = df1.index.str.lower()
	df2.index = df2.index.str.lower()
else:
	pass
if df1.empty == True:
    if TYPE == "i":
        result= pd.concat([df1, df2], axis=1, join='inner')
        result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+".txt", sep="\t", header=True)
    elif TYPE == "o":
        result = pd.concat([df1, df2], axis=1, join='outer') # join union
        #for traditional combining of both names, use the following
        result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+".txt", sep="\t", header=True)
    elif TYPE == "df1":
        result = pd.merge([df1, df2], how='left') # exact join from first dataframe
        result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+".txt", sep="\t", header=True)
    else:
        print ("Need TYPE {i= inner/intersection, o= outer/union, df1= exact join from first dataframe")
    
else:
    df1x= df1[df1.columns[0]]
    df1x = df1x.dropna(axis=0)
    df1y = df1[df1.columns[1:x]]
    df1= pd.concat([df1x, df1y], axis=1, join='inner') #join intersection after getting rid of NAs for class


    if DF3 != "none":
        df3 = pd.read_csv(DF3, sep='\t', index_col = 0)
        df3_name = os.path.basename(DF3)
        if LOWER == "T":
            df3.index.values = df3.index.values.str.lower()
        elif LOWER.lower() == "true":
            df3.index.values = df3.index.values.str.lower()
        else:
            pass
        col_list3= list(df3.columns.values)
        z= len(col_list3)
        if TYPE == "i":
            df2x = df2[df2.columns[1:y]]
            df3x = df3[df3.columns[1:z]]   
            df= pd.concat([df1, df2x, df3x], axis=1, join='inner') #join intersection    
            #added for combining binary and continuous data but same classes
            #df1name = str(DF1).split(".")[0]
            df.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+ str(df3_name)+".txt", sep="\t", header=True)
        elif TYPE == "o":
            result = pd.concat([df1, df2, df3], axis=1, join='outer') # join union
            #for traditional combining of both names, use the following
            result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+ str(df3_name)+".txt", sep="\t", header=True)
        elif TYPE == "df1":
            df2x = df2[df2.columns[1:y]]
            df3x = df3[df3.columns[1:z]]
            result = pd.concat([df1, df2x, df3x], axis=1, join_axes=[df1.index]) # exact join from first dataframe
            result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+ str(df3_name)+".txt", sep="\t", header=True)
        else:
            print ("Need TYPE {i= inner/intersection, o= outer/union, df1= exact join from first dataframe")

    else:
        if TYPE == "i":
            df2x = df2[df2.columns[1:y]]    
            df= pd.concat([df1, df2x], axis=1, join='inner') #join intersection    
            #added for combining binary and continuous data but same classes
            #df1name = str(DF1).split(".")[0]
            #df.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+".txt", sep="\t", header=True)
            df.to_csv(path_or_buf=str(df1_name)+"_ZF_TFs.txt", sep="\t", header=True)
        elif TYPE == "o":
            result = pd.concat([df1, df2], axis=1, join='outer') # join union
            #for traditional combining of both names, use the following
            result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+".txt", sep="\t", header=True)
        elif TYPE == "df1":
            result = pd.merge(df1, df2, how='left')  # exact join from first dataframe
            result.to_csv(path_or_buf=str(df1_name)+"_"+ str(df2_name)+".txt", sep="\t", header=True)
        else:
            print ("Need TYPE {i= inner/intersection, o= outer/union, df1= exact join from first dataframe")