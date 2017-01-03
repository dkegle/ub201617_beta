import glob
import itertools

import numpy
import pandas
import re
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram


class HierarchicalClustering:
    gene_distances = {}
    strains = []
    distances = None
    Z = None
    """
    An (n−1) by 4 matrix Z is returned. At the i-th iteration, clusters with indices Z[i,
    0] and Z[i, 1] are combined to form cluster n+i. A cluster with an index less than n
    corresponds to one of the n original observations. The distance between clusters Z[i,
    0] and Z[i, 1] is given by Z[i, 2]. The fourth value Z[i, 3] represents the number of
    original observations in the newly formed cluster."""

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

    def avg_linkage(self, c1, c2):
        return numpy.mean([self.distances.iloc[x[0], x[1]] for x in itertools.product(c1, c2)])

    def closest_clusters(self):
        return min([(self.avg_linkage(*map(self.indexes, c)), c) for c in
                    itertools.combinations(self.clusters, 2)], key=lambda x: x[0])

    def do_clustering(self):
        self.calculate_distances()
        i = 0
        while len(self.clusters) != 1:
            dist, (c1, c2) = self.closest_clusters()

            self.Z[i, 0] = self.cluster_names[str(c1)]
            self.Z[i, 1] = self.cluster_names[str(c2)]
            self.Z[i, 2] = dist
            self.Z[i, 3] = len(self.indexes([c1, c2]))

            print(c1, c2)
            self.clusters.remove(c2)
            self.clusters[self.clusters.index(c1)] = [c1 + c2]
            self.cluster_names[str([c1 + c2])] = next(self.cluster_idx_generator)
            i += 1

    def indexes(self, iterable):
        _l = [x for x in iterable if type(x) == int]
        for i in iterable:
            if hasattr(i, '__iter__'):
                _l.extend(self.indexes(i))
        return _l


hc = HierarchicalClustering()
hc.do_clustering()

countries = {}
for filename in glob.glob('data/genome/*.gb'):
    with open(filename, 'r') as f:
        country = re.search(r'/country="([a-zA-Z ,:]+)"', f.read()).group(1)
        countries[filename[12:-3]] = country

labels = [countries[str(hc.strains[x])] for x in range(len(hc.strains))]

plt.clf()
dendrogram(hc.Z, orientation='right', labels=labels)

label_colors = {'Guinea': 'r', 'Sierra Leone': 'k',
                'Gabon': 'b', 'Liberia': 'm', 'Democratic':'g'}

ax = plt.gca()
ylbls = ax.get_ymajorticklabels()
for lbl in ylbls:
    text = lbl.get_text()
    for key in label_colors:
        if key in text:
            lbl.set_color(label_colors[key])

plt.show()
