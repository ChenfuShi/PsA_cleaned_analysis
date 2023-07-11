import pandas as pd
import numpy as np

def quantile_normalize(df):
    """
    input: dataframe with numerical columns
    output: dataframe with quantile normalized values
    https://cmdlinetips.com/2020/06/computing-quantile-normalization-in-python/#:~:text=Quantile%20Normalization%20in%20Python%20When,in%20analyzing%20high%2Ddimensional%20datasets.
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)