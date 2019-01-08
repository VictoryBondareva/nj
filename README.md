# Построение филогенетического дерева методом ближайших соседей
Пример использования

```python
from nj import NJ_tree

names = ['a', 'b', 'c', 'd', 'e']
matrix = [[0], [5, 0], [9, 10, 0], [9, 10, 8, 0], [8, 9, 7, 3, 0]]
t = NJ_tree()
print(t.create_tree(names, matrix))
t.visualize()
```
