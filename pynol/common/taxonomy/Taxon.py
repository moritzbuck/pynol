class Taxon():

    def __repr__(self):
        repr_string = "{name} : {level} --- {string}"
        return repr_string.format(name = self.name, level = Taxon.rev_leves[self.level], string = self.get_tax_string())

    rev_leves = {
    0 : "root",
    1 : "domain",
    2 : "phylum",
    3 : "class",
    4 : "order",
    5 : "family",
    6 : "genus",
    7 : "species"
    }

    prefixes = {
    0 : "",
    1 : "d__",
    2 : "p__",
    3 : "c__",
    4 : "o__",
    5 : "f__",
    6 : "g__",
    7 : "s__"
    }

    levels = {
    "root" : 0,
    "domain" : 1,
    "phylum" : 2,
    "class" : 3,
    "order" : 4,
    "family" : 5,
    "genus" : 6,
    "species" : 7,
    "d" : 1,
    "p" : 2,
    "c" : 3,
    "o" : 4,
    "f" : 5,
    "g" : 6,
    "s" : 7
    }

    def __init__(self, level, name, parent):
        self.name = name
        self.parent = parent
        self.level = level if type(level) == int else Taxon.levels[level]

    def get_tax_string(self, full = False):
        if self.level > 1:
            return self.parent.get_tax_string(full) + ";" + (Taxon.prefixes[self.level] if full else "") + self.name
        else :
            return (Taxon.prefixes[self.level] if full else "") + self.name

    def is_child(self, parent):
        if self == parent:
            return True
        elif self.level == 0:
            return False
        else :
            return self.parent.is_child(parent)
