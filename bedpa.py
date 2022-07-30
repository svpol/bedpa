import pandas as pd
from pathlib import Path


class Bed6:

    def __init__(self, filepath):
        self.filepath = filepath
        self.directory = Path(filepath).parent
        self.stem = Path(filepath).stem.split('.')[0]
        self.df = pd.read_table(filepath, header=None, sep='\t')

    def get_chr_names(self):
        """
        Shows all chromosomes or scaffolds present in the object's data frame (and the initial file).
        :return: an array object with chromosome names.
        """
        return self.df[0].unique()

    def select_chr(self, chr_name):
        """
        Selects all records belonging to a given chromosome/scaffold and writes them to a bed file.
        :param chr_name: str, chromosome name, one of the values returned by get_chr_names.
        :return: path to the output file.
        """
        chr_selected = self.df[self.df[0] == chr_name]
        out_file = self.stem + "_" + chr_name + ".bed"
        out_path = Path(self.directory, out_file)
        chr_selected.to_csv(out_path, sep='\t', header=False, index=False)
        return out_path

    def select_strand(self, strand):
        """
        Selects all records for a specified DNA strand - plus or minus - and writes them to a bed file.
        :param strand: str, one of the values - 'plus' or 'minus'.
        :return: path to the output file.
        """
        if strand == "plus":
            selected_strand = self.df[self.df[5] == "+"]
        elif strand == "minus":
            selected_strand = self.df[self.df[5] == "-"]
        out_file = self.stem + "_" + strand + ".bed"
        out_path = Path(self.directory, out_file)
        selected_strand.to_csv(out_path, sep='\t', header=False, index=False)
        return out_path


class Bed12(Bed6):

    def __init__(self, filepath):
        super().__init__(filepath)

    def get_seq_type(self):
        """
        Shows all possible sequence types, e.g.: "exon", "gene", "UTR", etc.
        :return: an array object with all possible sequence types.
        """
        return self.df[7].unique()

    def count_seq_types(self):
        counts = self.df[7].value_counts()
        counts = counts.to_frame(name="count")
        share = self.df[7].value_counts() / len(self.df)
        share = share.to_frame(name="share")
        counts_share = pd.concat([counts, share], axis=1)
        counts_share = counts_share.rename_axis('seq_type').reset_index()
        return counts_share

    def select_seq_type(self, seq_type):
        """
        Selects all records of a particular type (like "exon", "gene", etc.) and writes them to a bed file.
        :param seq_type: str, sequence type, one of the values returned by get_seq_type.
        :return:
        """
        selected_types = self.df[self.df[7] == seq_type]
        out_file = self.stem + "_" + seq_type + ".bed"
        out_path = Path(self.directory, out_file)
        selected_types.to_csv(out_path, sep='\t', header=False, index=False)
        return out_path


if __name__ == "__main__":

    bed6 = Bed6('./samples/bed6/bed6.bed')
    print(bed6.get_chr_names())
    bed6_chr22_path = bed6.select_chr('chr22')
    del bed6
    bed6_chr22 = Bed6(bed6_chr22_path)
    bed6_chr22.select_strand("plus")
    del bed6_chr22

    bed12 = Bed12('./samples/bed12/bed12.bed')
    print(bed12.count_seq_types())
    print(bed12.get_chr_names())
    bed12_chr1_path = bed12.select_chr("chr1")
    del bed12
    bed12_chr1 = Bed12(bed12_chr1_path)
    print(bed12_chr1.get_seq_type())
    bed12_chr1.select_seq_type("exon")
    del bed12_chr1
