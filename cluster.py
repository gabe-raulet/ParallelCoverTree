import os
import sys
import numpy as np
from sklearn.neighbors import radius_neighbors_graph
from scipy.sparse.csgraph import connected_components

def read_vecs_file(fname):
    if fname.endswith(".fvecs"): dtype = np.float32
    elif fname.endswith(".ivecs"): dtype = np.int32
    elif fname.endswith(".bvecs"): dtype = np.uint8
    else: raise ValueError(f"Unknown file type: {fname}")

    filesize = os.path.getsize(fname)

    with open(fname, "rb") as f:
        d = int.from_bytes(f.read(4), byteorder="little")
        n = filesize // (4*(d+1))
        data = np.zeros((n,d), dtype=dtype)
        f.seek(0,0)
        for i in range(n):
            assert d == int.from_bytes(f.read(4), byteorder="little")
            data[i] = np.frombuffer(f.read(4*d), dtype=dtype, count=d)

    return data

def write_vecs_file(fname, data):
    if data.dtype == np.float32 and not fname.endswith(".fvecs"): fname += ".fvecs"
    elif data.dtype == np.int32 and not fname.endswith(".ivecs"): fname += ".ivecs"
    elif data.dtype == np.uint8 and not fname.endswith(".bvecs"): fname += ".bvecs"
    else: raise ValueError(f"Unknown vector type: {str(data.dtype)}")

    n, d = data.shape
    with open(fname, "wb") as f:
        for i in range(n):
            f.write(d.to_bytes(4, byteorder="little"))
            f.write(data[i].tobytes())

def distance(p, q):
    return np.sqrt(np.dot(p - q, p - q))

def get_clusters_from_labels(labels):
    n_comps = max(labels) + 1
    clusters = []
    for i in range(n_comps):
        clusters.append(np.where(labels==i)[0].tolist())
    return clusters

def get_clusters(data, epsilon):
    graph = radius_neighbors_graph(data, epsilon, mode="connectivity", include_self=True)
    n_comps, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
    return get_clusters_from_labels(labels)

def write_clusters(f, clusters):
    for cluster in clusters:
        f.write(" ".join((str(v) for v in cluster)) + "\n")

def main(points_fname, epsilon):
    data = read_vecs_file(points_fname)
    clusters = get_clusters(data, epsilon)
    write_clusters(sys.stdout, clusters)
    sys.stdout.flush()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <points> <epsilon>\n")
        sys.stderr.flush()
        sys.exit(1)
    main(sys.argv[1], float(sys.argv[2]))
