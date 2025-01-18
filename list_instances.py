import networkx as nx
import pandas as pd
import os

inputDir = "inputs/"

def max_grado(grafo):
    return max([grafo.degree(v) for v in grafo])

columns = ["name", "vertices", "degree"]
data = []
for input in os.listdir(inputDir):
    name = inputDir + input
    try:
        graph = nx.read_graphml(name)
    except Exception:
        continue
    n = graph.number_of_nodes()
    degree = max_grado(graph)
    data.append((input, n, degree))

df = pd.DataFrame(data, columns=columns)
df = df.sort_values(by="vertices", ignore_index=True)
df.to_csv('instances.csv')
