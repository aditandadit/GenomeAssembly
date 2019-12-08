# python2
import sys

"""
GENOME ASSEMBLY FROM K-MER
Task: Given the K-mer composition of some string,
      Perform Genome Assembly and return circular genome
      from which the k_mer came.

Genome Assembly => Shortest circular string that has all strings as substring
Input => Each of the 5396 lines of the input contains a single k-mer.
        The k-mer are given to you in alphabetical order because their true order is hidden from you.
        Each k-mer is 10 nucleotides long.
"""

'''
NOTES -
1. Started with Overlap Graph which uses Hamiltonian Cycle - No polynomial time solution exists
2. Now Using de bruijn Graphs - Eulerian Cycle - Linear Time solution exists
Hierholzer's Algorithm: O(Edges) Time
    Step 1: Check if Eulerian Cycle Exists => IN Degree = OUT Degree at each node
    Step 2: Choose any vertex V and follow a trail of edges till you return to V
            2 Cases
            Case 1 - All edges are Covered => DONE
            Case 2 - Some vertices on the tour have unexplored edges
    Step 3: Choose any vertex from previous tour that have unexplored edges
            and repeat this search process
            Join the previous and new tour
            
Implementation Details:
*   Maintain a Set of unused edges
*   Maintain a List of vertices on current tour that have unused edges

'''

'''
Genome Assembly =>
Step 1 : Build De brujin Graph from k-mers
Step 2 : Run Eulerian Cycle on this Graph to determine genome
'''


class EulerianCycle:
    def __init__(self, adj):
        self.n = len(adj)
        self.numUnexploredEdges = 0  # number of unexplored edges
        self.nodesWUE = dict()  # key: node with unused edges; value: the position of such node in the current path
        self.inDeg = dict()
        self.outDeg = dict()
        self.adjCurPos = dict()
        self.path = []
        self.unbalancedNode = []
        self.adjList = adj
        self.update_adj_list()

    def update_adj_list(self):
        for w, vList in self.adjList.items():
            self.inDeg[w] = self.inDeg.get(w, 0)
            for v in vList:
                self.inDeg[v] = self.inDeg.get(v, 0) + 1
            l = len(vList)
            self.outDeg[w] = l
            self.numUnexploredEdges += l
            self.adjCurPos[w] = 0

    def read_input(self):
        data = list(sys.stdin.read().strip().split())
        curMax = 0
        for i in range(len(data) // 3):
            curMax = max(int(data[i * 3]), curMax, max(list(map(int, data[i * 3 + 2].split(',')))))
        self.n = curMax + 1
        self.adjList = [[]] * self.n
        self.unusedEdges = [[]] * self.n
        self.inDeg = [0] * self.n
        self.outDeg = [0] * self.n
        self.adjCurPos = [0] * self.n
        for i in range(len(data) // 3):
            curIn = int(data[i * 3])
            self.adjList[curIn] = list(map(int, data[i * 3 + 2].split(',')))
            for v in self.adjList[curIn]:
                self.inDeg[v] += 1
            l = len(self.adjList[curIn])
            self.outDeg[curIn] = l
            self.numUnexploredEdges += l

    def add_edge(self):
        if type(self.adjList) is dict:
            for v in self.adjList.keys():
                if self.inDeg[v] != self.outDeg[v]:
                    if self.inDeg[v] < self.outDeg[v]:
                        self.unbalancedNode.append(v)
                    else:
                        self.unbalancedNode.insert(0, v)
            if len(self.unbalancedNode) > 0:
                self.adjList[self.unbalancedNode[0]].append(self.unbalancedNode[1])
                self.outDeg[self.unbalancedNode[0]] += 1
                self.inDeg[self.unbalancedNode[1]] += 1
            return
        for v in range(self.n):
            if self.inDeg[v] != self.outDeg[v]:
                if self.inDeg[v] < self.outDeg[v]:
                    self.unbalancedNode.append(v)
                else:
                    self.unbalancedNode.insert(0, v)
        if len(self.unbalancedNode) > 0:
            self.adjList[self.unbalancedNode[0]].append(self.unbalancedNode[1])
            self.outDeg[self.unbalancedNode[0]] += 1
            self.inDeg[self.unbalancedNode[1]] += 1
        return

    def explore(self, s):
        self.path.append(s)
        curr_pos = self.adjCurPos[s]
        cur_max_pos = self.outDeg[s]
        while curr_pos < cur_max_pos:
            self.adjCurPos[s] = curr_pos + 1
            if curr_pos + 1 < cur_max_pos:
                self.nodesWUE[s] = len(self.path) - 1
            else:
                if s in self.nodesWUE:
                    del self.nodesWUE[s]
            v = self.adjList[s][curr_pos]
            self.path.append(v)
            s = v
            curr_pos = self.adjCurPos[s]
            cur_max_pos = self.outDeg[s]
            self.numUnexploredEdges -= 1
        return

    def update_path(self, start_pos):
        l = len(self.path) - 1
        self.path = self.path[start_pos:l] + self.path[:start_pos]
        for node, pos in self.nodesWUE.items():
            if pos < start_pos:
                self.nodesWUE[node] = pos + l - start_pos
            else:
                self.nodesWUE[node] = pos - start_pos
        return

    def build_eulerian_cycle(self):
        if type(self.adjList) is dict:
            w, vList = self.adjList.popitem()
            self.adjList[w] = vList
            self.explore(w)
        else:
            self.explore(0)

        while self.numUnexploredEdges > 0:
            node, pos = self.nodesWUE.popitem()
            self.update_path(pos)
            self.explore(node)
        return self.path

    def print_path(self):
        print('->'.join([str(node) for node in self.path]))


class GenomeAssembly_k_mer_composition:
    def __init__(self):
        self.k, self.adj = self.read_data()
        self.path = EulerianCycle(self.adj).build_eulerian_cycle()
        print(self.reconstruct_from_path(self.path)[:-self.k + 1])

    def read_data(self):
        data = list(sys.stdin.read().strip().split())
        adj = self.de_brujin(len(data[0]), data)
        return len(data[0]), adj

    @staticmethod
    def de_brujin(k, patterns):
        adjdb = dict()
        for p in patterns:
            if p[:k - 1] in adjdb:
                adjdb[p[:k - 1]].append(p[1:])
            else:
                adjdb[p[:k - 1]] = []
                adjdb[p[:k - 1]].append(p[1:])
            if p[1:] not in adjdb:
                adjdb[p[1:]] = []
        return adjdb

    @staticmethod
    def reconstruct_from_path(path):
        return path[0] + ''.join(seq[-1] for seq in path[1:])


if __name__ == "__main__":
    GenomeAssembly_k_mer_composition()

'''
What value of K to use? => SPADES Program uses K = 31 to 127
Genome Bases are not randomly spread out, so K-mers of large size repeat frequently
'''
