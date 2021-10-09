# curl -L 'https://api.genome.ucsc.edu/getData/track?genome=hg19;chrom=chr8;start=128745680;end=128755674;track=hub_124277_RLFS'
# curl -L 'https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr8;start=128745680;end=128755674;'

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

rloopr_verts, rloopr_positions = load_bedgraph("./myc/chr8_hub_124277/d5194783-f902-4282-9579-5cb895706be2.bedgraph", 0.35)


# QMRLFS
qmrlfs = GenomeBrowserRecord(0.15)
qmrlfs.load_file("./myc/chr8_hub_124277/RLFS.txt")

# hub_124277_DRIP1_peaks_NT2_Sanz_2016
#drip1 = GenomeBrowserRecord(0.3)
#drip1.load_file("./myc/chr8_hub_124277/hub_124277_DRIP1_peaks_NT2_Sanz_2016.txt")

# hub_124277_Fibroblast_Nadel_2015
#rdip = GenomeBrowserRecord(0.4)
#rdip.load_file("./myc/chr8_hub_124277/hub_124277_Fibroblast_Nadel_2015.txt")


# chr11_neat1_hub_124277_RNA-seq_r2_NT2_Sanz_2016
#nt2 = GenomeBrowserRecord(0.9)
#nt2.load_file("./neat1/chr11_neat1_hub_124hub_124277_Fibroblast_Nadel_2015
print(qmrlfs.chromosome)
print(qmrlfs.track_name)


fig, ax = plt.subplots()

patch = patches.PathPatch(Path(qmrlfs.verts, codes * int(len(qmrlfs.verts)/2)), edgecolor='red', lw=4)
ax.add_patch(patch)
ax.text(qmrlfs.track_start + 20, 0.2, "QmRLFS-finder", color='red')

patch = patches.PathPatch(Path(rloopr_verts, codes * int(len(rloopr_verts)/2)), edgecolor='green', lw=4)
ax.add_patch(patch)
ax.text(qmrlfs.track_start + 2, 0.4, "Rloop tracker"+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, rloopr_positions) +" %)", color='green')

#patch = patches.PathPatch(Path(drip1.verts, codes * int(len(drip1.verts)/2)), edgecolor='black', lw=2)
#ax.add_patch(patch)
#ax.text(drip1.track_start, 0.26, drip1.track_name+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, drip1.positions) +"%)", color='black')

#atch = patches.PathPatch(Path(rdip.verts, codes * int(len(rdip.verts)/2)), edgecolor='black', lw=2)
#ax.add_patch(patch)
#ax.text(rdip.track_start, 0.36, rdip.track_name+" (similarity to Qm-RFLS results: "+ cmp_similarity(qmrlfs.positions, rdip.positions) +"%)", color='black')



ax.set_xlim(qmrlfs.track_start, qmrlfs.track_end)
ax.set_ylim(0, 0.6)

plt.xlabel("Sequence position [bp]")
plt.ylabel("R-loop analysis type")
plt.show()
