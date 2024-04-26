import os
import sys
import numpy as np
from snnpy import *
import time

dim = 2

def read_points(fname):
    t = -time.perf_counter()
    with open(fname, "rb") as f:
        n = int.from_bytes(f.read(8), byteorder="little")
        data = np.zeros((n,dim), dtype=np.float32)
        for i in range(n):
            data[i] = np.frombuffer(f.read(4*dim), dtype=np.float32, count=dim)
    t += time.perf_counter()
    print(f"[time={t:.4f}] :: (read_points) [filename='{fname}']")
    return data

def build_graph(data, radius):
    t = -time.perf_counter()
    snn_model = build_snn_model(data)
    t += time.perf_counter()
    print(f"[time={t:.4f}] :: (build_snn_model) [n={data.shape[0]}]")
    snn_time = t

    n_edges = 0
    edges = []
    t = 0
    for i in range(len(data)):
        t1 = -time.perf_counter()
        ind = snn_model.query_radius(data[i], radius, return_distance=False)
        t1 += time.perf_counter()
        t += t1
        n_edges += len(ind)
        edges.append(ind)
    print(f"[time={t:.4f}] :: (queries)")
    query_time = t
    print(f"[time={query_time + snn_time:.4f}] :: (build_graph) [n_verts={data.shape[0]},n_edges={n_edges},avg_deg={n_edges/data.shape[0]:.4f},radius={radius:.3f}]")
    return edges, n_edges

def write_graph(edges, n_edges, ofname):
    t = -time.perf_counter()
    with open(ofname, "w") as f:
        f.write(f"{len(edges)} {n_edges}\n")
        for u in range(len(edges)):
            for v in sorted(edges[u]):
                f.write(f"{u+1} {v+1}\n")
    t += time.perf_counter()
    print(f"[time={t:.4f}] :: (write_graph) [filename='{ofname}']")

def main(ifname, radius, ofname=None):
    data = read_points(ifname)
    edges, n_edges = build_graph(data, radius)
    if ofname: write_graph(edges, n_edges, ofname)
    return 0

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <points> <radius> [graph]\n")
        sys.stderr.flush()
        sys.exit(1)

    ofname = None
    if len(sys.argv) >= 4: ofname = sys.argv[3]

    sys.exit(main(sys.argv[1], float(sys.argv[2]), ofname))
