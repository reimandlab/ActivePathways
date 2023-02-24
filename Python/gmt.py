import numpy as np

class GMT:
    def __init__(self, name_dict, genes_dict):
        self.name_dict = name_dict
        self.genes_dict = genes_dict

def read_gmt(infile):
    # read file and create a dictionary
    gmt = [line.strip("\n") for line in open(infile) if len(line.strip("\n")) > 2]
    
    ids = list(map(lambda x: x.split("\t")[0], gmt))
    names = list(map(lambda x: x.split("\t")[1], gmt))
    genes = list(map(lambda x: np.array(x.split("\t")[2:]), gmt))

    # dictionary of ids : names / gene arrays
    name_dict = dict(zip(ids, names))
    genes_dict = dict(zip(ids, genes))

    return GMT(name_dict, genes_dict)

def write_gmt(gmt, filename):
    # write gmt to file
    with open(filename, "w") as outfile:
        for key in gmt.name_dict.keys():
            outfile.write("\t".join([key, gmt.name_dict[key], '\t'.join(gmt.genes_dict[key])]) + "\n")

def is_gmt(gmt):
    return isinstance(gmt,GMT)

def make_background(gmt):
    # return a unique list of genes from the gmt
    try:
        genes = []
        [genes.extend(gene_list) for gene_list in gmt.genes_dict.values()]
        return np.unique(genes)

    except:
        print("gmt is not a valid GMT object")
        return
