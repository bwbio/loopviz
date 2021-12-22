# -*- coding: utf-8 -*-
"""
@author: bw98jh
"""

import pandas as pd
import numpy as np
import math
from tqdm import tqdm

#return df with [chr, centroid1, centroid2, metric]
def preprocess(path, metric='fdrDonut'):
    df = pd.read_csv(path,sep='\t')[1:]
    df = df.rename(columns = {'#chr1':'chr'})    
    df = df[df.chr == df.chr2].drop(columns='chr2')
    df.chr = df.apply(lambda x: 'chr'+str(x.chr), axis=1)
    df = df[['chr','centroid1','centroid2',metric]]
    df[['centroid1','centroid2']] = df[['centroid1','centroid2']].astype(int)
    return df

#return dataframes df_list = [df_unique0, df_static, df_unique1]
def differential(f0,f1,metric='fdrDonut',mindist=10000): 
    
    df0, df1 = preprocess(f0), preprocess(f1)
    
    k = 10**7
    
    chroms = sorted(df0.chr.unique())
    
    static = []
    unique0,unique1 = [],[]
    
    for chrom in chroms:
        df0chr, df1chr = df0[df0.chr == chrom], df1[df1.chr == chrom]
        maxbp = max(df0chr.centroid2.max(), df1chr.centroid2.max())
        roundbins =  int(math.ceil(maxbp/k))
        
        overhang = (df0chr.centroid2-df0chr.centroid1).max() + 500000
        
        bins = [((i-1)*k - overhang, i*k + overhang) for i in range(roundbins)]
        
        
        dropper0, dropper1 = [],[] #collect common loop indices
        
        for L,R in tqdm(bins,leave=False, desc=chrom):
            df0seg, df1seg = df0chr[(df0chr.centroid1>L) & (df0chr.centroid2<R)],\
                             df1chr[(df1chr.centroid1>L) & (df1chr.centroid2<R)]
        
            for i, row0 in df0seg.iterrows():
                c1, c2 = row0.centroid1, row0.centroid2
                matching = df1seg[((abs(df1seg.centroid1-c1) <= mindist)) &\
                                   (abs(df1seg.centroid2-c2) <= mindist)]
                
                if len(matching) > 0:
                    dropper0.append(i)
                    dropper1 = dropper1+matching.index.to_list()
        
        dropper0, dropper1 = list(set(dropper0)), list(set(dropper1))
        
        unique0.append(df0chr.drop(index=dropper0)) #exclude common loop indices
        unique1.append(df1chr.drop(index=dropper1))
        static.append(df0chr.loc[dropper0])
        
    df_unique0 = pd.concat(unique0)
    df_unique1 = pd.concat(unique1)
    df_static = pd.concat(static)
     
    df_list = [df_unique0, df_static, df_unique1]
    return df_list

#convert df to .interact format output dataframe
def convert_df(df, samplename, hex_color, metric):
    
    out = pd.DataFrame(columns = ['#chrom','chromStart','chromEnd',
                                  'name','score','value','exp','color',
                                  'sourceChrom','sourceStart','sourceEnd',
                                  'sourceName','sourceStrand',
                                  'targetChrom','targetStart','targetEnd',
                                  'targetName','targetStrand'
                                  ])
    
    for i in ['#chrom','sourceChrom','targetChrom']: out[i] = df.chr
    for i in ['chromStart','sourceStart','sourceEnd']: out[i] = df.centroid1
    for i in ['chromEnd','targetStart','targetEnd']: out[i] = df.centroid2
    
    for i in ['sourceName','sourceStrand','targetName','targetStrand']: out[i] = '.'
    out.name = samplename+out.index.astype(str)
    out.exp = samplename
    out.score = 0
    out.value = df[metric]
    out.color = hex_color
    return out

#wrapper for convert_df
def collate_dfs(df_list, metric):
    label_list = ['lost','static','gained']
    color_list = ['#ff003c','#000000','#42b6f5']
    out_list = []
    
    for i in range(len(df_list)):
        label = label_list[i]
        color = color_list[i]
        df    = df_list[i]
        
        out_list.append(convert_df(df,label,color,metric))
    
    out = pd.concat(out_list)
    return out

#write output from collate_dfs to .interact file
def write_output(out, samplename, outpath):
    header = ' '.join(['track type=interact',
                       'track name="{}"'.format(samplename),
                       'name="{}"'.format(samplename),
                       'description="{}"'.format(samplename),
                       'interactDirectional=false',
                       'maxHeightPixels=200:100:50',
                       'visibility=full',
                       ]
                      )+'\n' 
    
    out.to_csv(outpath, sep='\t', index=False)
    with open(r'{}'.format(outpath), 'r+') as f: 
      base = f.read() 
      with open(r'{}'.format(outpath), 'w+') as f: 
          f.write(header + base)   

#full pipeline
def loopviz(f0,f1,samplename,outpath='loopviz_test.interact',metric='fdrDonut',mindist=10000):
    if outpath.endswith('.interact'):
        df_list = differential(f0,f1,metric,mindist)
        out = collate_dfs(df_list, metric)
        write_output(out,samplename,outpath)
        return out
    else:
        print('Failed: please end [outpath] with ".interact"')

#%% Testing

#%% Testing

# f0 = 'HK1_VC_merged_loops.bedpe'
# f1 = 'HK1_TH_merged_loops.bedpe'
# metric = 'fdrDonut'

# x=loopviz(f0,f1,samplename='HK1_THZ1',outpath='HK1_THZ1_data.interact')

#%% Values

# string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records" 
# uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region" 
# uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
# string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
# uint score;          "Score (0-1000)"
# double value;        "Strength of interaction or other data value. Typically basis for score"
# string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
# string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
# string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
# uint sourceStart;    "Start position in chromosome of source/lower/this region"
# uint sourceEnd;      "End position in chromosome of source/lower/this region"
# string sourceName;   "Identifier of source/lower/this region"
# string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
# string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
# uint targetStart;    "Start position in chromosome of target/upper/this region"
# uint targetEnd;      "End position in chromosome of target/upper/this region"
# string targetName;   "Identifier of target/upper/this region"
# string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"