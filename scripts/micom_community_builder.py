import micom
from micom import Community
from os import makedirs
from os.path import isfile
from micom.workflows import workflow
import pandas as pd
from micom import load_pickle
from sys import argv

import numpy as np

logger = micom.logger.logger
try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 30

makedirs("models", exist_ok = True)

'''
enter the file to create pickle files
'''
taxonomy = pd.read_csv(argv[1]).query("relative > 1e-3")

taxonomy["file"] = taxonomy.file.apply(
    lambda ids: ["../../arrivale/agora_models/" + i for i in ids.split("|")]
)
taxonomy["name"] = taxonomy.genus
assert not taxonomy.name.str.contains(" ").any()
taxonomy = taxonomy.rename(columns={"name" : "id", "reads": "abundance"})

diet = pd.read_csv("../../micom/western_diet.csv")
diet.index = diet.reaction = diet.reaction.str.replace("_e", "_m")
diet = diet.flux * diet.dilution


def build_and_save(args):
    try:
        s, tax,prefix = args
        filename = "models/" + prefix+ "_" + s + ".pickle"
        if isfile(filename):
            return
        com = Community(tax, id=s, progress=False)
        ex_ids = [r.id for r in com.exchanges]
        logger.info(
            "%d/%d import reactions found in model.",
            diet.index.isin(ex_ids).sum(),
            len(diet),
        )
        com.medium = diet[diet.index.isin(ex_ids)]
        com.to_pickle(filename)
    except:
        print('problem with', s)
        
samples = taxonomy.samples.unique()
args = [(s, taxonomy[taxonomy.samples == s],argv[2]) for s in samples]
workflow(build_and_save, args, max_procs)
