import glob
import itertools

import pandas
import re

import matplotlib.pyplot as plt
import networkx, pylab
from networkx.drawing.nx_agraph import graphviz_layout
from Bio import Phylo
from Bio.Phylo import TreeConstruction


def max_value(inputlist):
    return max([max(sublist) for sublist in inputlist])

#What color to give to the edges?
e_color = '#ccccff'
#What colors to give to the nodes with similar labels?
color_scheme = {'Guinea': 'r', 'Sierra Leone': 'k',
                'Gabon': 'y', 'Liberia': 'm', 'Democratic':'g', 'Nigeria':'b'}
#What sizes to give to the nodes with similar labels?
size_scheme = {'PLACEHOLDER':200}

#Edit this to produce a custom label to color mapping
def label_colors(label):
	color_to_set = 'blue'
	for label_subname in color_scheme:
		if label_subname in label:
			color_to_set = color_scheme[label_subname]
	return color_to_set

#Edit this to produce a custom label to size mapping
def label_sizes(label):
	#Default size
	size_to_set = 20
	for label_subname in size_scheme:
		if label_subname in label:
			size_to_set = size_scheme[label_subname]
	return size_to_set

    
    
class NeighbourJoining:
    gene_distances = {}
    strains = []
    distances = None
    tree = None
    countries = None

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

    def calculate_distances(self):
        def strain_distances(g1, g2):
            # Sum of all gene distances for 2 strains
            return sum(self.gene_distances[gene][g1, g2] for gene in self.gene_distances.keys())

        self.distances = pandas.DataFrame(index=self.strains, columns=self.strains)
        for x, y in itertools.permutations(range(len(self.strains)), 2):
            self.distances.ix[x, y] = strain_distances(self.strains[x], self.strains[y])        
    
        
    def neigh_join(self):
        # Fix the given DataFrame
        self.calculate_distances()
        self.distances.fillna(0, inplace=True)

        tmp = [x[:i+1] for i, x in enumerate(self.distances.values.tolist())]
        max_val = max_value(tmp)
        norm = list(map(lambda x: [float(i)/max_val for i in x], tmp))

        dm = TreeConstruction._DistanceMatrix(names=self.distances.keys().tolist(), matrix=norm)
        
        constructor = TreeConstruction.DistanceTreeConstructor()
        tree = constructor.nj(dm)
        tree.ladderize()
        self.tree = tree


    def plot(self):
        G = Phylo.to_networkx(self.tree)

        node_sizes = []
        labels = {}
        node_colors = []

        countries = {}
        for filename in glob.glob('genes/*.gb'):
            with open(filename, 'r') as f:
                country = re.search(r'/country="([a-zA-Z ,:]+)"', f.read()).group(1)
                countries[filename[6:-3]] = country.rsplit(':')[0]
        countries_names = [countries[str(self.strains[x])] for x in range(len(self.strains))]

        for n in G:
            label = str(n)
            if 'Inner' in label:
                # These are the inner tree nodes -- leave them blank and with very small sizes.
                node_sizes.append(1)
                labels[n] = ''
                node_colors.append(e_color)
            else:
                label = countries[str(n)]
                # Size of the node depends on the labels!
                node_sizes.append(140)
                # Set colors depending on our color scheme and label names
                node_colors.append(label_colors(label))
                # set the label that will appear in each node
                # labels[n] = label

        pos = graphviz_layout(G)

        f = plt.figure(1,figsize=(12, 12))
        ax = f.add_subplot(1, 1, 1)
        for label in set(countries_names):
            ax.plot([0], [0],
                    color=label_colors(label),
                    label=label)
        plt.legend(scatterpoints=1, loc='upper right', shadow=False, prop={'size': 8})

        networkx.draw_networkx(G, pos, edge_color=e_color, node_size=node_sizes, labels=labels, with_labels=True,
                               node_color=node_colors, ax=ax)

        pylab.savefig('neighbour_joining.png')
        pylab.show()



if __name__ == "__main__":
    nj = NeighbourJoining()
    nj.neigh_join()
    nj.plot()
