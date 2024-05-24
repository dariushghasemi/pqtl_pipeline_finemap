from pathlib import Path
import pandas as pd


# Define input for the rules
data = []
with open(config["path_sumstats"], "r") as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    seqid = ".".join(p.stem.split(".")[:3])
    data.append((seqid, str(p)))

analytes = (
    pd.DataFrame.from_records(data, columns=["seqid", "path_sumstats"])
    .set_index("seqid", drop=False)
    .sort_index()
)


def get_sumstats(wildcards):
    return analytes.loc[wildcards.seqid, "path_sumstats"]


# define the functions generating files' path
def ws_path(file_path):
    return str(Path(config.get("path_workspace"), file_path))

# function to take the meta-analysis loci
def sc_path(file_path):
    return str(Path(config.get("path_loci"), file_path))