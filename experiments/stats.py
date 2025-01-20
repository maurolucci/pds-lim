import pandas as pd
import os

def get_stats(statsPath: str) -> pd.DataFrame:
    columns = ["solver", "instance", "vertices", "edges", "propagating_vertices", "omega",
              "variables", "constraints", "run", "lower_bound", "upper_bound", "gap",
              "result", "nodes", "t_solver", "callbacks", "t_callback", "lazy_constraints"]
    data = []
    for fileName in os.listdir(statsPath):
        with open(statsPath + "/" + fileName) as file:
            for line in file:
                s = line.split(",")
                data.append((s[0], s[1], int(s[2]), int(s[3]), int(s[4]), int(s[5]),
                             int(s[6]), int(s[7]), int(s[8]), float(s[9]), int(s[10]),
                             float(s[11]), s[12], int(s[13]), int(s[14]), int(s[15]),
                             int(s[16]), int(s[17])))
    return pd.DataFrame(data, columns)

def main(statsPath: str, statsName: str) -> None:
    df = get_stats(statsPath)
    df.to_csv(statsName)
