import networkx as nx
from revolutionhtl.nhxx_tools import get_nhx


class TreePolytomies:
    def __init__(self, og: int, tree: nx.DiGraph, X: dict[int, dict[int, list[int | str]]]):
        self._og = og
        self._tree = tree
        self._X = X.copy()
        # X = {
        #       x: {
        #           y_1: C_1 = [z_1, z_2,...]
        #           y_2: C_2 = [z_1, z_2,...]
        #           ...
        #       },
        #
        #       x': {
        #           y'_1: C'_1 = [z'_1, z'_2,...]
        #           y'_2: C'_2 = [z'_1, z'_2,...]
        #           ...
        #       },
        #       ...
        # }

    def get_og(self) -> int:
        return self._og

    def get_tree(self) -> nx.DiGraph:
        return self._tree

    def get_nodes_with_polytomies(self) -> list[int]:
        return list(self._X.keys())

    def get_ys(self, x: int) -> list[int]:
        return list(self._X[x].keys())

    def get_cluster(self, x, y_i) -> list[int | str]:
        return self._X[x][y_i]

    def __str__(self):
        text = f"OG = {self.get_og()}\n" \
               f"\ttree = {self.get_tree()}\n" \
               f"\tnodes with polytomies = {self.get_nodes_with_polytomies()}\n"

        for x in self.get_nodes_with_polytomies():
            text += f"\t{x = }\n"
            for y_i in self.get_ys(x):
                text += f"\t\t{y_i = }: C_i = {self.get_cluster(x, y_i)}\n"

        return f"{text}\tNewick: {get_nhx(self.get_tree(), name_attr='label')}\n"
