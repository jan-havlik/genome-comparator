import json

seq = open("Neat1_chr8:128745680-128755674_raw", "r")
data = json.load(seq)

with open("Neat1_chr8:128745680-128755674", "w") as out:
    out.write(data["dna"].upper())
