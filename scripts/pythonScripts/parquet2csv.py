import argparse
import os
import re
import pyarrow.parquet as pq
import pandas as pd

## functions

def replace_with_dict(value,chrMap):
    """Parquet file has weird chrom codes so 
    need to replace w/ numbers"""
    if value in chrMap:
        return chrMap[value]
    else:
        return None

def processPq(inFile,chrMap,dataDir):
    """Takes in pq file, changes chrom code to ID and filters out low quality (failed)
    reads w/ card < 2"""
    df = pq.read_table(source=f'{dataDir}/{inFile}').to_pandas()
    filtered_df = df[df['pass_filter'] == True]
    filtered_df.loc[:,'chrom'] = filtered_df['chrom'].apply(replace_with_dict,chrMap = chrMap)
    filtered_df = filtered_df.dropna(subset=['chrom'])
    filtered_df.reset_index(drop=True, inplace=True)
    readCounts = filtered_df.groupby('read_idx').size().reset_index(name='count')
    filteredRead_ix = readCounts[readCounts['count'] >= 2]
    filtered_df = filtered_df[filtered_df['read_idx'].isin(filteredRead_ix['read_idx'])]
    return(filtered_df)

def writeOutFile(inFile,chrMap,dataDir,outDir):
    outName = re.sub("_unphased.at.pore_c.parquet",".csv.gz",inFile).replace("run","").replace("_batch","")
    print(outName)
    outDF = processPq(inFile,chrMap,dataDir)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outDF.to_csv(f'{outDir}/{outName}',index = False,compression = 'gzip')

def main():
    ## Setup
    print("Processing input arguments")
    args = parse_args()
    rootDir = args.root
    dataDir = f"{rootDir}/{args.dataDir}"
    outDir = f"{rootDir}/{args.outDir}"
    chromMap = f"{rootDir}/{args.chromMap}"
    chrMap = pd.read_table(chromMap).set_index('alias2')['chrom'].to_dict()

    print("Processing input files")
    allFiles = os.listdir(dataDir)
    for inFile in allFiles:
        writeOutFile(inFile,chrMap,dataDir,outDir)
    

def parse_args():
    parser = argparse.ArgumentParser(description="""Takes in a parquet file and converts
                                     to a csv format amenable to previously written
                                     bash preprocessing script""")

    parser.add_argument("root",type=str, help="root directory")
    parser.add_argument("dataDir",type=str, help="input data dir relative to root")
    parser.add_argument("outDir",type=str, help="out dir relative to root")
    parser.add_argument("chromMap",type=str, help="chromosome map relative to root")

    return parser.parse_args()

if __name__ == "__main__":
    main()