import pandas as pd
import numpy as np
import sqlite3

def mapping_func(grouped_peptide_df):
    grouped_peptide_df = grouped_peptide_df[["ID", "MODIFIED_SEQUENCE"]]
    grouped_peptide_df["ID_TYPE"] = grouped_peptide_df["MODIFIED_SEQUENCE"].apply(lambda x: "UNIMOD_ID" if "UniMod" in x else "CODENAME_ID")
    grouped_peptide_df = grouped_peptide_df[["ID_TYPE", "ID"]]
    
    if len(grouped_peptide_df) == 1:
        if grouped_peptide_df["ID_TYPE"].item() == "CODENAME_ID":
            grouped_peptide_df = pd.concat([grouped_peptide_df, pd.DataFrame({"ID": [-1], "ID_TYPE":"UNIMOD_ID"})])
        else:
            grouped_peptide_df = pd.concat([grouped_peptide_df, pd.DataFrame({"ID": [-1], "ID_TYPE":"CODENAME_ID"})])
    
    grouped_peptide_df = pd.pivot_table(grouped_peptide_df, columns="ID_TYPE")
    return grouped_peptide_df

def make_mapping_table(osw_file):
    peptide_df = pd.read_sql_query("SELECT * FROM PEPTIDE", osw_file)
    peptide_df["tmp"] = peptide_df["MODIFIED_SEQUENCE"].replace("\\(\\w+\\)|\\(\\w+[:]\\d+\\)", "(@)", regex=True)
    peptide_df = peptide_df.groupby("tmp")
    peptide_df = test_df.apply(mapping_func) 
    peptide_df.columns.name = None
    peptide_df = peptide_df.reset_index(drop=True)

    return peptide_df

# peptide_df.to_sql("UNIMOD_CODENAME_MAPPING", osw_file,index=False) 
# write above line to write output of make_mapping_table to the osw file