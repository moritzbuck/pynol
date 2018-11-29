from pynol.common.taxonomy.Taxon import Taxon
from tqdm import tqdm


class Taxonomy():

    def __getitem__(self, key):
        if key[-4:] == ";s__" :
            key = key.replace(";s__", "")
        first = self.taxa_dict.get(key)
        if first:
            return first
        else :
            taxa = key.split(";")
            assert all([self.taxa_dict.get(t) for t in taxa]), "This is either not a taxon string or some taxa don't exist in this DB"
            assert self.check_consistency(key), "This is not a consistent taxon string"
            return self.taxa_dict.get(taxa[-1])

    def __init__(self):
        self.taxa_dict = {}
        self.root = Taxon(0,"root", None)
        self.taxa_dict["root"] = self.root

    def from_gtdb_taxonomy_file(self, file):
        with open(file) as handle:
            all_data  = {l.strip().split("\t")[0] : l.strip().split("\t")[1] for l in handle}

        all_strings = set(all_data.values())

        for s in tqdm(all_strings):
            parent = self.root
            fields = [(f.split("__")[0], f.split("__")[1]) for f in s.split(";")]
            for f in fields :
                key = f[0] + "__" + f[1]
                if not (f[0] == 's' and f[1] == ''):
                    if self.taxa_dict.get(key):
                        parent = self.taxa_dict.get(key)
                    else :
                        self.taxa_dict[key] = Taxon(f[0], f[1], parent)
                        parent = self.taxa_dict[key]
        self.gtdb = { k : self[v] for k,v in all_data.items()}


    def check_consistency(self , tax_string):
        parent = self.root
        taxa = tax_string.split(";")
        if taxa[-1] == "s__":
            del taxa[-1]
        for t in taxa:
            child = self[t]
            if not child.is_child(parent):
                return False
            parent = child
        return True

    def check_monophyly(self, taxa_list):
        tax_set ={ i : set([self[Taxon.prefixes[i] + t.get_parent(i)] for t in taxa_list if t.get_parent(i) ]) for i in range(7)}
        if len({ k : v for k,v in tax_set.items() if len(v) >1}) > 0:
            return False
        else :
            return True
