import glob
import itertools

import numpy
import pandas

class CreateDistanceMatrix:
    gene_distances = {}
    strains = []
    distances = None
    Z = None

    def __init__(self):
        for filename in glob.glob('distance_matrix-*.txt'):
            with open(filename, 'r') as f:
                gene = filename[16:-4]
                self.gene_distances[gene] = {}
                for line in f.readlines():
                    if line.startswith('#'):
                        continue
                    g1, g2, distance = line.split(' ')
                    self.strains.append(g1)
                    self.gene_distances[gene][g1, g2] = int(distance.strip('\n'))
                    self.gene_distances[gene][g2, g1] = int(distance.strip('\n'))

        self.strains = list(set(self.strains))
        self.clusters = [[x] for x in range(len(self.strains))]
        self.cluster_idx_generator = iter(range(len(self.strains), 123123123123))
        self.cluster_names = {str(self.clusters[i]): i for i in range(len(self.clusters))}
        self.Z = numpy.ndarray((len(self.strains) - 1, 4), dtype=numpy.double)

    def calculate_distances(self):
        def strain_distances(g1, g2):
            # Sum of all gene distances for 2 strains
            return sum(self.gene_distances[gene][g1, g2] for gene in self.gene_distances.keys())

        self.distances = pandas.DataFrame(index=self.strains, columns=self.strains)
        for x, y in itertools.permutations(range(len(self.strains)), 2):
            self.distances.ix[x, y] = strain_distances(self.strains[x], self.strains[y])
        self.distances.fillna(0, inplace=True)

    def to_csv(self, filename):
        self.calculate_distances()
        with open(filename, 'w') as f:
            self.distances.to_csv(f)
        


if __name__ == "__main__":
    dm = CreateDistanceMatrix()
    dm.to_csv("distance_matrix.csv")
