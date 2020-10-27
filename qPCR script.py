import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

def graphqpcr(csv):
    """Takes qPCR csv output file with one reference gene and three targets and returns graphs"""

    #Open qPCR data and remove unneeded data
    df = pd.read_csv(csv, header = 7)
    df = df.drop(df.index[-4:-1])
    df = df.drop(df.columns[0:1], axis = 1)
    df = df.drop(df.columns[2:8], axis = 1)
    df = df.drop(df.columns[5:], axis = 1)
    df = df.dropna()
    df = df.set_index('Sample Name')
    
    #Create reference gene dataframe
    ref = (df.loc[df['Target Name'] == df.iloc[0,0]])
    
    #Create testing genes dataframe
    gene1 = (df.loc[df['Target Name'] == df.iloc[3,0]])
    gene2 = (df.loc[df['Target Name'] == df.iloc[6,0]])
    gene3 = (df.loc[df['Target Name'] == df.iloc[9,0]])
    genes = gene1.append(gene2)
    genes = genes.append(gene3)
    
    #Isolate CT mean values
    genes_CT = genes['Cт Mean']
    ref_CT = ref['Cт Mean']
    ref_CT = ref_CT.append(ref['Cт Mean'])
    ref_CT = ref_CT.append(ref['Cт Mean'])

    #Calculate delta CT values
    delta_mean = genes_CT - ref_CT
    calcs = 2**(-delta_mean)

    #Remove duplicates to make the values easier to play with
    trim_calcs = calcs.drop_duplicates(keep = 'first')

    #Isolate standard deviations
    genes_SD = genes['Cт SD']
    ref_SD = ref['Cт SD']
    ref_SD = ref_SD.append(ref['Cт SD'])
    ref_SD = ref_SD.append(ref['Cт SD'])

    #Calculate SE values as per analysis spreadsheet
    delta_mean_SD = np.sqrt(genes_SD**2 + ref_SD**2)
    delta_SE = delta_mean_SD/np.sqrt(3)
    CT_max = delta_mean + delta_SE
    CT_min = delta_mean - delta_SE
    calcs_max = 2**(-CT_max)
    calcs_min = 2**(-CT_min)
    SE = (calcs_min - calcs_max)/2

    #Remove duplicates
    trim_SE = SE.drop_duplicates(keep = 'first')
  

    #Define the values for each graph
    cat = [df.index[0], df.index[12], df.index[24], df.index[36]]
    val1 = trim_calcs[0:4]
    val2 = trim_calcs[4:8]
    val3 = trim_calcs[8:12]
    SE1 = trim_SE[0:4]
    SE2 = trim_SE[4:8]
    SE3 = trim_SE[8:12]

    #Set figure size
    %pylab inline
    plt.figure(figsize(30,10))
    plt.subplots_adjust(wspace = 0.5)

    #Make the compound graphs
    plt.subplot(1,3,1)
    plt.bar(cat, val1, yerr = SE1, capsize = 5, color = 'grey', edgecolor = 'black')
    plt.xticks(fontsize = 30)
    plt.yticks(fontsize = 25)
    plt.ylabel(df.iloc[3,0]+'/'+df.iloc[0,0], fontsize = '45')

    plt.subplot(1,3,2)
    plt.bar(cat, val2, yerr = SE2, capsize = 5, color = 'grey', edgecolor = 'black')
    plt.xticks(fontsize = 30)
    plt.yticks(fontsize = 25)
    plt.ylabel(df.iloc[6,0]+'/'+df.iloc[0,0], fontsize = '45')

    plt.subplot(1,3,3)
    plt.bar(cat, val3, yerr = SE3, capsize = 5, color = 'grey', edgecolor = 'black')
    plt.xticks(fontsize = 30)
    plt.yticks(fontsize = 25)
    plt.ylabel(df.iloc[9,0]+'/'+df.iloc[0,0], fontsize = '45')
   
    plt.savefig(df.index[12] + ':' + df.iloc[3,0] + ',' + df.iloc[6,0] + ',' + df.iloc[9,0] + '.png')


    plt.show()
