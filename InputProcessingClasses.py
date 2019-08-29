import pandas as pd
import os
import subprocess as sb
from subprocess import Popen
import requests
import sys


class GwasSignal:

    """
    This class transform UK data.
    Input: path to UKBB format file.
    Output:
    """

    def __init__(self, path, chrom=0, MAF=0.05):
        self.path = path
        self.chrom = chrom
        self.tmp = os.path.join('tmp', 'gwas_signal.tsv')
        self.MAF = MAF
        self.pval = 10**-8

    @staticmethod
    def convert_37_to_38(df):
        """
        This requires the most of the time.
        """
        server = "https://rest.ensembl.org"
        for index in df.index:
            rs = df['MarkerName'].loc[index]
            ext = f"/variation/human/{rs}?"
            r = requests.get(
                server + ext, headers={"Content-Type": "application/json"})
            if not r.ok:
                # r.raise_for_status()
                print('Error: SNP coordinates can not be converted')
                continue
            decoded = r.json()
            location = decoded['mappings'][0]['start']
            df['Position_GRCh38'].loc[index] = location

        return df

    @staticmethod
    def cut_variants_with_little_MAF(df, MAF):
        filtered_df = df[((df['AC']) / (df['nCompleteSamples'] * 2)) >= MAF]
        filtered_df.index = range(len(filtered_df.index))
        return filtered_df

    @staticmethod
    def cut_variants_with_high_pval(df, pval):
        filtered_df = df[df['pval'] <= pval]
        filtered_df.index = range(len(filtered_df.index))
        return filtered_df

    def filter_df(self, df, MAF, pval):
        return self.cut_variants_with_high_pval(
            self.cut_variants_with_little_MAF(df, MAF), pval)

    def change_file_structure(self, path):
        df = self.filter_df(pd.read_csv(path, sep='\t'), self.MAF, self.pval)

        new_df = pd.DataFrame(
            {'Chromosome': pd.Series(list(map(lambda x: x.split(':')[0], df['variant']))),
             'Position_GRCh38': pd.Series(list(map(lambda x: x.split(':')[1], df['variant']))),
             'MarkerName': df['rsid'],
             'Effect_allele': pd.Series(list(map(lambda x: x.split(':')[3], df['variant']))),
             'Non_Effect_allele': pd.Series(list(map(lambda x: x.split(':')[2], df['variant']))),
             'Beta': df['beta'],
             'SE': df['se'],
             'Pvalue': df['pval']})
        return new_df

    def cut_chrom(self, df):
        if self.chrom != 0:
            df = df[df['Chromosome'] == str(self.chrom)]
        return df

    def make_new_df(self):

        new_df = self.convert_37_to_38(
            self.cut_chrom(
                self.change_file_structure(self.path)))

        new_df.to_csv(self.tmp, sep='\t', index=False)
        # os.remove(self.tmp)
        return os.path.abspath(self.tmp)


class PostGap:

    """
    This class launch postgap and write output to:
    1) stderr - errors during postgap activation
    2) stdout - prints during postgap
    3) postgap_output
    """

    def __init__(self, path, disease, postgap_tmp=True):
        self.path = path
        self.disease = disease
        self.db = os.path.join('databases_dir', 'databases')
        self.tmp = os.path.join(os.path.split(self.path)[0],
                                'postgap_output.tsv')
        self.postgap_tmp = postgap_tmp

    @staticmethod
    def write_postgap_tmp(pop_obj, path):
        stdout, stderr = pop_obj.communicate()
        with open(os.path.join(path, 'stdout.txt'), 'w') as out:
            with open(os.path.join(path, 'stderr.txt'), 'w') as err:
                print(stdout.decode("utf-8"), file=out)
                print(stderr.decode("utf-8"), file=err)

    def get_postgap_data(self):
        os.chdir('postgap-master/')
        try:
            process = Popen('python2 POSTGAP.py --summary_stats {0} --disease {1} --database_dir {2} --output {3}'.format(
                self.path, self.disease, self.db, self.tmp), shell=True, stdout=sb.PIPE, stderr=sb.PIPE)
            os.chdir("..")
            if self.postgap_tmp == True:
                self.write_postgap_tmp(process, 'tmp')
            return pd.read_csv(self.tmp, sep='\t')
        except FileNotFoundError:
            print('Postgap Error! See stderr file.')

    def remove_tmp(self):
        os.remove(self.tmp)
