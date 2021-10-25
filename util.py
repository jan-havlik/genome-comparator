import json

seq = open("chr11_raw.txt", "r")
data = json.load(seq)

with open("Neat1_chr11:65188245-65215011", "w") as out:
    out.write(data["dna"].upper())
