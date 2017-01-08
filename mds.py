import glob
import itertools
import re
import numpy
import pandas
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn import preprocessing

class Mds:
    gene_distances = {}
    strains = []
    distances = None
    pos = None

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

    def normalize_data(self):
        self.distances = ( self.distances- self.distances.mean()) / (self.distances.max() - self.distances.min())

    def do_mds(self):
        self.calculate_distances()

        #min_max_scaler = preprocessing.MinMaxScaler()

        self.distances.fillna(0, inplace=True)
        #normalized = min_max_scaler.fit_transform(self.distances)
        #self.distances = pandas.DataFrame(normalized)


        print(self.distances)
        clf = manifold.MDS(n_components=2, max_iter=3000, eps=1e-6, dissimilarity="precomputed", n_jobs=1)
        self.pos = clf.fit(self.distances).embedding_

        # Rotate the data
        # clf = PCA(n_components=2)

        # self.pos = clf.fit_transform(self.distances)

mds = Mds()
mds.do_mds()

countries = {}
for filename in glob.glob('data/genome/*.gb'):
    with open(filename, 'r') as f:
        country = re.search(r'/country="([a-zA-Z ,:]+)"', f.read()).group(1)
        countries[filename[12:-3]] = country

countrys = [countries[str(mds.strains[x])] for x in range(len(mds.strains))]

label_colors = {'Guinea': 'r', 'Sierra Leone': 'gold',
                'Gabon': 'c', 'Liberia': 'orange', 'Democratic Republic of the Congo':'b',
                'Nigeria':'g'}
# plt.scatter(X_true[:, 0], X_true[:, 1], color='navy', s=s, lw=0,
#             label='True Position')

plotData = dict()
for label, pos in zip(mds.strains, mds.pos):
    country = countries[label]
    key = country[:country.find(':')] if country.find(':') > 0 else country
    if key in plotData:
        plotData[key][0].append(pos[0])
        plotData[key][1].append(pos[1])
    else:
        plotData[key] = ([pos[0]],[pos[1]])

for key in plotData:
    print(plotData[key])
    plt.scatter(plotData[key][0], plotData[key][1], color=label_colors[key], s=120, lw=0, label=key, alpha=0.7)

plt.legend(scatterpoints=1, loc='upper left', shadow=False, prop={'size':8})




plt.show()