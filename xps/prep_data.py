import numpy
import csv
import re
import glob

import pdb

def get_pos_k(k, counts):
    offsets = {"F": {5: 0, 10: 3.5, 20: 9},
               "C": {5: 15, 10: 18.5, 20: 24},
               "B": {5: 30, 10: 33.5, 20: 37}}
    pos = 0
    if k[0] == "E":
        if k[1] == 5:
            pos = -4
        elif k[1] == 10:
            pos = -3
        else:
            pos = -2
    else:
        othk = (k[0], k[1], k[2], "a" if k[3]=="s" else "s")
        if othk in counts:
            pos = 2.5*(k[2]=="D") + (k[3]=="a") + offsets[k[0]][k[1]]
        else:
            pos = 1.5*(k[2]=="D") + offsets[k[0]][k[1]]
    return pos

        
counts = {}
for f in glob.glob("data_simul_*.csv"):
    num_lines = sum(1 for line in open(f))
    tmp = re.match("data_simul_(?P<series>(?P<group>[A-Z])(?P<sparsity>[DS])(?P<direction>[sa])(?P<poolsize>[0-9]+))(?P<polar>[OPN])\.csv", f)
    k = (tmp.group("group"), int(tmp.group("poolsize")), tmp.group("sparsity"), tmp.group("direction"))
    if k not in counts:
        counts[k] = {}
    counts[k][tmp.group("polar")] = num_lines

pos = dict([(k, get_pos_k(k, counts)) for k in counts.keys()])

print("%%% xticks")
print(sorted(pos.values()))
print("%%% series nodes")
print("\n".join(["\\node (n%s%s%s%s) at (axis cs:%s,0) {};%%" %(k[0], k[2], k[3], k[1], pos[k]) for k,c in sorted(counts.items(), key=lambda x: pos[x[0]])]))
print("%%% labels loop")
print("\n".join(["{n%s%s%s%s}/{%s}/{%s}/{%s}/{%s}/{%s}/{%s}/{%s}/{%s},%%" %(k[0], k[2], k[3], k[1], pos[k], k[0], k[1], k[2], k[3], c.get("O", "-"), c.get("P", "-"), c.get("N", "-")) for k,c in sorted(counts.items(), key=lambda x: pos[x[0]])]))

# data from diagnosis_cut.csv plotting "cooccur" (top), "mu" (middle) and "fitfilt" (bottom) grouped by "run_name" (columns) and "sim" (colors)
#### Simul recall boxplots
file_in="diagnosis_cut.csv"
file_out="data_simul-%s.csv"
vls = ("cooccur", "mu", "fitfilt")
dt_map = {}
with open(file_in) as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in csvreader:
        #### polarity of associations
        if float(row["sim"].strip("'")) == 1:
            assoc = "P"
        elif float(row["sim"].strip("'")) == -1:
            assoc = "N"
        else:
            assoc = "O"
        #### direction: symmetric/assymmetric
        direction = "a" if row["run_name"][-1] == "+" else "s"
        #### poolsize, nb of species
        if re.search("Sp5", row["run_name"]) is not None:
            poolsize = "5"
        elif re.search("Sp10", row["run_name"]) is not None:
            poolsize = "10"
        else:
            poolsize = "20"
        #### group, env only, facilitation, competition, both
        if re.search("FacComp", row["run_name"]) is not None:
            group = "B"
        elif re.search("Fac", row["run_name"]) is not None:
            group = "F"
        elif re.search("Comp", row["run_name"]) is not None:
            group = "C"
        else:
            group = "E"
        #### density: dense/sparse
        density = "S" if re.search("Sparse", row["run_name"]) is not None else "D"            
        k = "".join([group,density,direction,poolsize,assoc])
        if k not in dt_map:
            dt_map[k] = []
        dt_map[k].append(tuple([float(row[v].strip("'")) for v in vls]))
# print([(k,len(v)) for (k,v) in dt_map.items()])

for (k,vs) in dt_map.items():    
    with open(file_out % k, "w") as fo:
        fo.write("\n".join([",".join(["%f" % vx for vx in v]) for v in vs]))


#### Simul recall scatter plot
file_in="assoc_group.csv"
file_out="data_simulStats.csv"
map_vs = {"group": {"F": "F", "C": "C", "FC": "B"},
          "assoc": {"1": "P", "'-1": "N"},
          "direction": {"1": "a", "0": "s"},
          "density": {"Sparse": "S", "Dense": "D"}}
fields = [("group", "%s"), ("density", "%s"), ("direction", "%s"), ("assoc", "%s"),
          ("poolsize", "%d"), ("support", "%d"),
          ("precision", "%f"), ("recall", "%f"), ("f1-score", "%f")]
vrs, fmt = zip(*fields)
vrs = list(vrs)
vrs[-1] = "fone"
fo = open(file_out, "w")
fo.write(",".join(vrs)+"\n")
with open(file_in) as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in csvreader:
        vout = []
        for (v,f) in fields:
            if v in map_vs:
                tmp = map_vs[v].get(row[v].strip(), "?")
            else:
                tmp = float(row[v].strip("'"))
                if v == "recall":
                    tmp /= 100
            vout.append(f % tmp)
        fo.write(",".join(vout)+"\n")
fo.close()


#### MI histograms
file_in="raw_AlpsMI.txt"
file_out="data_AlpsMI.csv"
nb_vars = 8
A = numpy.loadtxt(file_in)
X = numpy.reshape(A, (-1,nb_vars))[3:-2]
Y = numpy.concatenate((1*[numpy.arange(X.shape[1])+1], X))
numpy.savetxt(file_out, Y.T, fmt='%f', delimiter=",")

#### Alps HSM weights boxplots
file_in="weights.csv"
file_out="data_AlpsHSMw.csv"
var = "variable"
val = "value"
dt_map = {}
with open(file_in) as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in csvreader:
        if row[var] not in dt_map:
            dt_map[row[var]] = []
        dt_map[row[var]].append(float(row[val].strip("'")))
hs = sorted(dt_map.keys())
Y = numpy.array([dt_map[h] for h in hs])
numpy.savetxt(file_out, Y.T, fmt='%f', delimiter=",", header=",".join(hs))


#### Alps plants scores
file_in="scores.csv"
file_out="data_AlpsPlantScoresTMP.csv"
fo = open(file_out, "w")
vrs = ["species", "genus","score", "prevalence", "nboccur", "avgcountnz", "avgcount"]
fo.write(",".join(vrs+["gcl", "pngl","slbl"])+"\n")
with open(file_in) as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in csvreader:
        Gchr = chr(ord("A")+(ord(row["genus"][0])-ord("A"))%11)
        fo.write(",".join([row[v] for v in vrs]+[Gchr,"90",row["species"][0:2]+row["species"][4:7]])+"\n")
fo.close()

#### Alps plants scores
file_in="scores.csv"
file_out="data_AlpsPlantScoresTMP.csv"
fo = open(file_out, "w")
vrs = ["species", "genus","score", "prevalence", "nboccur", "avgcountnz", "avgcount"]
fo.write(",".join(vrs+["gcl", "pngl"])+"\n")
with open(file_in) as csvfile:
    csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in csvreader:
        Gchr = chr(ord("A")+(ord(row["genus"][0])-ord("A"))%8)
        fo.write(",".join([row[v] for v in vrs]+[Gchr,"90"])+"\n")
fo.close()
