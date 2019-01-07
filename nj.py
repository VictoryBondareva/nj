import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo import BaseTree 
from Bio.Phylo import Consensus
import copy
import itertools

class NJ_tree:
	def __init__(self):
		self.tree = None

	def vizualize(self):
		if self.tree != None:
			Phylo.draw(self.tree)

	def create_tree(self, names, matrix):
		if not names or not matrix:
			return self.tree

		distance_matrix = DistanceMatrix(names, matrix)
		dm = copy.deepcopy(distance_matrix)
		clades = [BaseTree.Clade(None, name) for name in dm.names]

		# print(clades)
		# print(dm)
		# clades[0].clades.append(clades[2])
		while len(clades) != 1:
			q_matrix = []
			q_names = dm.names
			for ind_i, i in enumerate(dm.matrix):
				tmp = []
				for ind_j, j in enumerate(i):
					if ind_i == ind_j:
						tmp.append(0)
						continue
					tmp.append((len(dm) - 2)*j - sum(dm[ind_i]) - sum(dm[ind_j]))
				q_matrix.append(tmp)
			q_matrix = DistanceMatrix(q_names, q_matrix)
			# print(q_matrix)

			min_i = float('Inf')
			min_j = float('Inf')
			q_min = float('Inf')
			for ind_i, i in enumerate(q_matrix):
				for ind_j, j in enumerate(i):
					if j < q_min:
						q_min = j
						min_i = ind_i
						min_j = ind_j
			# print(q_min, min_i, min_j)

			if len(clades) == 2:
				# print('c:', clade_j)
				if min_i == 0:
					clade_j = clades[min_j]
					clade_j.branch_length = dm[min_i][min_j]
					clades[min_i].clades.append(clade_j)
					del clades[min_j]
					# print(clades)
					break
				if min_i == 1:
					clade_i = clades[min_i]
					clade_i.branch_length = dm[min_i][min_j]
					clades[min_j].clades.append(clade_i)
					del clades[min_i]
					# print(clades)
					break


			dist_i = 0.5*dm[min_i][min_j] + (sum(dm[min_i]) - sum(dm[min_j]))/(2*(len(dm) - 2))
			dist_j = dm[min_i][min_j] - dist_i
			# print(dist_i, dist_j)

			clade_i = clades[min_i]
			clade_j = clades[min_j]
			clade_i.branch_length = dist_i
			clade_j.branch_length = dist_j

			tmp_clade = BaseTree.Clade(None, dm.names[min_i] + dm.names[min_j])
			tmp_clade.clades.append(clade_i)
			tmp_clade.clades.append(clade_j)

			clades[min_j] = tmp_clade
			del clades[min_i]
			# print(clades)

			tmp_dist = []
			for k in range(len(dm)):
				if k == min_j or k == 0:
					tmp_dist.append(0)
					continue
				tmp_dist.append(0.5*(dm[min_i][k] + dm[min_j][k] - dm[min_i][min_j]))
			dm[min_j] = tmp_dist

			dm.names[min_j] = dm.names[min_i] + dm.names[min_j]
			del dm[min_i]

		self.tree = BaseTree.Tree(clades[0], rooted = False)
		# print('res: ', BaseTree.Tree(clades[0], rooted = False))
		return self.tree



# constructor = DistanceTreeConstructor()
# names = ['a', 'b', 'c', 'd', 'e']
# matrix = [[0], [52, 0], [58, 13, 0], [90, 66, 11, 0], [83, 77, 2, 59, 0]]
# tree = constructor.nj(DistanceMatrix(None))
# print(tree)
# # Phylo.draw(tree)
		
# t = NJ_tree()
# print(t.create_tree(None))
# t.vizualize()

# if(Consensus._equal_topology(tree, t.tree)):
# 	print('Success!')


# for size in [3, 5, 7]:
# 	names = [str(i) for i in range(size)]
# 	matrix = []

# 	for i in range(size):
# 		tmp = [np.random.randint(0, 100) for j in range(i)]
# 		tmp.append(0)
# 		matrix.append(tmp)
	
# 	dm = DistanceMatrix(names, matrix)
# 	print(dm)
# 	tree = constructor.nj(dm)
# 	t = NJ_tree().create_tree(dm)


# 	print(tree)
# 	print(t)
# 	print(Consensus._equal_topology(tree, t))

# 	# pylab.figure(1)
# 	Phylo.draw(tree)
# 	# pylab.figure(2)
# 	Phylo.draw(t)
# 	# pylab.show()


