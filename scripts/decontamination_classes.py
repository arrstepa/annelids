"""
Classes for TPM-based cross-contamination filtration

This module allows user to:
- load TPM tables;
- calculate TPM assymetry between transcripts;
- visualize log2(TPMratio) distribution;
- identify probably contaminant transcripts.
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class ContigCounter():
    def __init__(self, contig_pairs: pd.DataFrame):
        self.contig_pairs = contig_pairs.copy()

    def __len__(self):
        return len(self.contig_pairs)
    
    def __repr__(self):
        return f"{self.__class__.__name__}(n_pairs={len(self)})"


class TPMProcessor():
    """
    TPM files processor: loads TPM tables and aggregates expression values with one of supported methods:
    - max
    - mean
    - median
    """

    def __init__(self, filepaths, sep='\t', aggregation_function = 'max'):
        if isinstance(filepaths, str):
            self.filepaths = [filepaths]
        else:
            self.filepaths = filepaths
        
        self.sep = sep
        self.data = None
        self.aggregation_function = aggregation_function
    
    
    def __repr__(self):
        return f"TPMProcessor(filepath={self.filepaths}, aggregation_function={self.aggregation_function})"

    def __call__(self):
        return self.load()

    def _aggregate(self, df):

        numeric_df = df.select_dtypes(include='number')

        if numeric_df.empty:
            raise ValueError('No numeric columns, invalid TPM file')


        if self.aggregation_function == 'max':
            return numeric_df.max(axis=1)
        elif self.aggregation_function == 'mean':
            return numeric_df.mean(axis=1)
        elif self.aggregation_function == 'median':
            return numeric_df.median(axis=1)
        else:
            raise ValueError('Unknown aggregation function. Try max, mean or median instead.')


    def load(self):
        dfs = []
        for path in self.filepaths:
            df = pd.read_csv(path, sep=self.sep)
            df = df.reset_index().rename(columns={'index': 'contigs'})
            df['tpm_value'] = self._aggregate(df)
            dfs.append(df[['contigs', 'tpm_value']])
        self.data = pd.concat(dfs, ignore_index=True)
        return self.data
    


class CrossContaminationDetector(ContigCounter):
    """
    Calculates log2(TPM ratio) between transcript pairs

    Large TPM asymetry is considered as contamination
    """
    def __init__(self, contig_pairs: pd.DataFrame, tpm_processor: TPMProcessor):
        super().__init__(contig_pairs)
        
        self.tpm_file = tpm_processor()
        self.contigs_tpm_dict = self.tpm_file.set_index('contigs')['tpm_value']

    
    def __call__(self):
        return self.count_log2ratio()
    

    def count_log2ratio(self):
        self.contig_pairs["tpm1"] = self.contig_pairs.iloc[:, 0].map(self.contigs_tpm_dict)
        self.contig_pairs["tpm2"] = self.contig_pairs.iloc[:, 1].map(self.contigs_tpm_dict)

        self.contig_pairs['log2_ratio'] = np.log2((self.contig_pairs['tpm1'] + 1 )/ 
                                                  (self.contig_pairs['tpm2'] + 1))
        return self.contig_pairs
        


class ThresholdFiltrator(ContigCounter):
    """
    Returns recomended threshold and visualisation
    """

    def __call__(self):
        return self.visualise()
        

    def visualise(self):
        plt.hist(np.abs(self.contig_pairs['log2_ratio'].dropna()), bins=100)
        plt.xlabel("|log2(TPM1 / TPM2)|")
        plt.ylabel("Count")
        plt.title("Absolute log2_ratio")
        plt.show()


    def recommend_threshold(self, cutoff=0.5):
        n = len(self.contig_pairs)
        list_of_thresholds = []

        for treshold in np.arange(0, 3, 0.5):
            contamination_frac = len(self.contig_pairs[np.abs(self.contig_pairs['log2_ratio']) >= treshold]) / n
            list_of_thresholds.append((treshold, contamination_frac))
       
        recommendation = min(
            (i for i in list_of_thresholds if i[1] > cutoff),
            key= lambda x: x[0],
            default=None
        )
        
        result_df = pd.DataFrame(list_of_thresholds, columns=['threshold', 'cross-contamination_frac'])
        return recommendation, result_df
    

    def filter(self, threshold):
        to_eliminate = pd.concat([self.contig_pairs[self.contig_pairs['log2_ratio'] >= threshold].iloc[:, 1],
                         self.contig_pairs[self.contig_pairs['log2_ratio'] <= -threshold].iloc[:, 0]],
                         ignore_index=True)
        return to_eliminate
