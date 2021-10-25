# curl -L 'https://api.genome.ucsc.edu/getData/track?genome=hg19;chrom=chr11;start=65188245;end=65194245;track=hub_124277_RLFS'
# curl -L 'https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr11;start=65188245;end=65194245;'

import json
from itertools import product, islice

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib import patches

VERTS = []
codes = [ # same code repeated multiple times (for each point)
    Path.MOVETO,
    Path.LINETO,
]

def load_bedgraph(path, pos):
    rloopr_verts = []
    rloopr_positions = []
    with open(path, 'r') as bed_f:
        for line in islice(bed_f, 2, None):  # skip first 2 lines of header
            split = line.split()
            # list of coordinates of every rloop
            rloopr_verts.append((float(split[1]), pos))
            rloopr_verts.append((float(split[2]), pos))
            rloopr_positions += range(int(split[1]), int(split[2]))
    return (rloopr_verts, rloopr_positions)


def cmp_similarity(ref, dst):
    total = len(set(ref+dst))
    return "%d" % (int(len(list(set(ref) & set(dst))) / total * 100))

class Rloop:

    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.model = ""
        self.strand = ""

    def __eq__(self, other):
        print(f"{self.start}-{self.end}\t{other.start}-{other.end}")
        if other.start in [self.start, self.start+1, self.start-1]:  # range of +- 1 to both sides
            if other.end in [self.end, self.end+1, self.end-1]:
                return True
        return False

class GenomeBrowserRecord:

    def __init__(self, graph_position):
        self.track_name = ""
        self.result_count = 0
        self.genome = ""
        self.chromosome = ""
        self.positions = []
        self.track_start = 0
        self.track_end = 0
        self.graph_position = graph_position
        self.verts = []

    def load_file(self, path):
        json_file = open (path, "r")
        data = json.loads(json_file.read())
        
        self.genome = data["genome"]
        self.track_name = data["track"]
        self.chromosome = data["chrom"]
        self.result_count = data["itemsReturned"]
        self.track_start = data["start"]
        self.track_end = data["end"]

        
        src = data[self.track_name] if isinstance(data[self.track_name], list) else data[self.track_name].get(self.chromosome)
        startString = "chromStart" if isinstance(data[self.track_name], list) else "start"
        endString = "chromEnd" if isinstance(data[self.track_name], list) else "end"

        for res in src:
            self.verts.append(([float(res[startString]), self.graph_position]))
            self.verts.append(([float(res[endString]), self.graph_position]))
            self.positions += range(int(res[startString]), int(res[endString]))


# QMRLFS
qmrlfs = GenomeBrowserRecord(0.05)
qmrlfs.load_file("./qmrlfs.txt")

# Rloop tracker
rloopr_verts, rloopr_positions = load_bedgraph("./chr11_rloopt.bedgraph", 0.15)

#G4 Hunter
g4_verts, g4_positions = load_bedgraph("./chr11_g4.bedgraph", 0.25)

# hub_124277_Fibroblast_rdip_Nadel_2015_peaks
rdip = GenomeBrowserRecord(0.35)
rdip.load_file("./Fibroblast_rdip_Nadel_2015_peaks.txt")
rdip.track_name = rdip.track_name.replace("hub_124277_", "")

# hub_124277_epithelial_rdip_Nadel_2015_peaks
rdip2 = GenomeBrowserRecord(0.45)
rdip2.load_file("./Epithelial_rdip_Nadel_2015_peaks.txt")
rdip2.track_name = rdip2.track_name.replace("hub_124277_", "")

print(qmrlfs.chromosome)
print(qmrlfs.track_name)


fig, ax = plt.subplots()

patch = patches.PathPatch(Path(qmrlfs.verts, codes * int(len(qmrlfs.verts)/2)), edgecolor='red', lw=4)
ax.add_patch(patch)
ax.text(qmrlfs.track_start + 20, 0.06, "QmRLFS-finder", color='red')

patch = patches.PathPatch(Path(rloopr_verts, codes * int(len(rloopr_verts)/2)), edgecolor='green', lw=4)
ax.add_patch(patch)
ax.text(qmrlfs.track_start + 2, 0.16, "Rloop tracker"+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, rloopr_positions) +" %)", color='green')

patch = patches.PathPatch(Path(g4_verts, codes * int(len(g4_verts)/2)), edgecolor='black', lw=2)
ax.add_patch(patch)
ax.text(qmrlfs.track_start + 2, 0.26, "G4 hunter"+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, g4_positions) +" %)", color='black')

patch = patches.PathPatch(Path(rdip.verts, codes * int(len(rdip.verts)/2)), edgecolor='black', lw=2)
ax.add_patch(patch)
ax.text(rdip.track_start, 0.36, rdip.track_name+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, rdip.positions) +"%)", color='black')

patch = patches.PathPatch(Path(rdip2.verts, codes * int(len(rdip2.verts)/2)), edgecolor='black', lw=2)
ax.add_patch(patch)
ax.text(rdip2.track_start, 0.46, rdip2.track_name+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, rdip2.positions) +"%)", color='black')

ax.set_xlim(qmrlfs.track_start, qmrlfs.track_end)
ax.set_ylim(0, 0.6)

plt.xlabel("Sequence position [bp]")
plt.ylabel("R-loop analysis type")
plt.show()
