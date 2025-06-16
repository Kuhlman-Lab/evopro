from typing import Any, List, Tuple, Dict
import copy
import networkx as nx

class Graph(nx.Graph):
    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f'<{list(self.nodes)}, {list(self.edges)}, {list(self.nodes[i]["v"] for i in self.nodes)}>'

    def bridges(self) -> List[Tuple[int, int]]:
        return list(nx.bridges(self))

    def relabel_nodes(self, rmap: Dict[int, int]) -> "Graph":
        return nx.relabel_nodes(self, rmap)

    def clear_cache(self) -> None:
        self._Data_cache = None


class GraphLike:
    # The Graph object whose attributes can be accessed.
    _graph: Graph

    def __init__(self, graph: Graph) -> None:
        self.graph = graph

    @property
    def graph(self) -> Graph:
        return self._graph

    @graph.setter
    def graph(self, new_graph: Graph) -> None:
        assert type(new_graph) is Graph
        self._graph = new_graph

    @property
    def nodes(self) -> List[Dict[str, Any]]:
        return self.graph.nodes

    @property
    def edges(self) -> List[Dict[str, Any]]:
        return self.graph.edges

    def __len__(self) -> int:
        return len(self._graph)

    def copy(self) -> "GraphLike":
        return copy.deepcopy(self)

    def remove_node(self, *args, **kwargs) -> None:
        self.graph.remove_node(*args, **kwargs)

    def remove_edge(self, *args, **kwargs) -> None:
        self.graph.remove_edge(*args, **kwargs)

    def reset_graph(self) -> None:
        self.graph = Graph()

    def initial_state(self) -> "GraphLike":
        """Return the intial state of the GraphLike. The initial state is a
        copy of the GraphLike with the underlying Graph reset.

        Returns:
            GraphLike: the initial state.
        """
        graphlike_copy = self.copy()
        graphlike_copy.reset_graph()
        return graphlike_copy