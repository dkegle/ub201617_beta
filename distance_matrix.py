from Bio import SeqIO
from Bio.Seq import Seq
from os import listdir
from itertools import combinations
from multiprocessing import Pool

genes = ["VP24", "NP", "VP35", "VP40", "VP30", "GP", "L"]
path = "xdata\chD_viruses_\data\genome\\"

def printGeneNamesAndSizes(seqobj):
   for feat in seqobj.features:
      if feat.type == "gene":
         print(" ".join([feat.qualifiers['gene'][0], "from", str(feat.location.start), "to", str(feat.location.end)]))

def getGene(seqobj, gene_name):
   for feature in seqobj.features:
      if feature.type == "gene" and feature.qualifiers['gene'][0] == gene_name:
         return seqobj[feature.location.start : feature.location.end]
   return Seq("")

def editDistance(s1, s2):
   def scoring_function(a, b):
      return 0 if a == b else 1

   # initialization
   M = {(0,0): 0}
   M.update({(0, j+1): j+1 for j in range(len(s1))})
   M.update({(i+1, 0): i+1 for i in range(len(s2))})

   # DP table
   for j in range(1, len(s1) + 1):
      for i in range(1, len(s2) + 1):
         score_up = M[i-1, j] + scoring_function(s2[i - 1], "-")
         score_diag = M[i-1, j-1] + scoring_function(s2[i - 1], s1[j - 1])
         score_left = M[i, j-1]+ scoring_function("-", s1[j - 1])
         M[i, j] = min(score_up, score_diag, score_left)

   return M[len(s2), len(s1)]

def getDistance(sample_1, sample_2, _genes):
   dist = 0
   for gene in _genes:
      gene_1 = getGene(sample_1["seqobj"], gene)
      gene_2 = getGene(sample_2["seqobj"], gene)
      dist += editDistance(gene_1, gene_2)
   return " ".join([sample_1["id"], sample_2["id"], str(dist)])

if __name__ == "__main__":
   samples = [{"seqobj": SeqIO.read(path + sample, "genbank"), "id": sample.rstrip(".gb")} for sample in listdir(path)]

   # testing
   #prvi = SeqIO.read(path + "436409269.gb", "genbank")
   #drugi = SeqIO.read(path + "824038961.gb", "genbank")
   #samples = [{"seqobj": prvi, "id": "436409269"}, {"seqobj": drugi, "id": "824038961"}]

   gene_subset = [genes[2]]   # only one gene atm
   out = open("distance_matrix.txt", "w")
   out.write("# Columns are as follows: sample_1 sample_2 distance")
   process_pool = Pool(processes = 8)
   results = [process_pool.apply_async(getDistance, args=(sample_1, sample_2, gene_subset)) 
         for (sample_1, sample_2) in combinations(samples, 2)]
   
   for result in results:
      out.write("\n" + result.get())
   out.close()

   # testing
   #printGeneNamesAndSizes(prvi)
   #gene_1 = getGene(prvi, "VP24")
   #gene_2 = getGene(drugi, "VP24")
   #print(editDistance(gene_1, gene_2))

   print("done")