import pandas as pd

def write_job_array_csv(path, df):
    df.to_csv(path)

def dict_from_row(row):
    return dict(zip(row.keys(), row))       

def rows_to_df(rs):
    ks = rs[0].keys()
    df = pd.DataFrame(columns=ks)
    for r in rs:
        df = pd.concat([df, pd.DataFrame(
            columns=ks,
            data=[dict_from_row(r).values()]
            )], ignore_index=True)
    return df
