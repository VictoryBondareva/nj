# Построение филогенетического дерева методом ближайших соседей
Пример использования

```python
from nj import NJ_tree
from Bio.Phylo.TreeConstruction import DistanceMatrix

names = ['a', 'b', 'c', 'd', 'e']
matrix = [[0], [52, 0], [58, 13, 0], [90, 66, 11, 0], [83, 77, 2, 59, 0]]
t = NJ_tree()
print(t.create_tree(names, matrix))
t.vizualize()
```
