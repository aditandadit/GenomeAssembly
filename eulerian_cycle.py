# python2
import sys

class EulerianCycle:
    def __init__(self):
        self.numExploredEdges = 0  # number of explored edges
        self.nodesWUE = dict()  # key: node with unused edges; value: the position of such node in the current path
        self.path = [] # final Path
        self.adjList = [] # Graph as Adjacency list

        # Check if Eulerian Path can Exist
        is_balanced = self.read_input()
        if is_balanced:
            print('1')
            self.build_eulerian_cycle()
        else:
            print('0')

        self.print_path()

    def explore(self, s):
        self.path.append(s)
        curPos = self.adjCurPos[s]
        curMaxPos = self.outDegree[s]
        while curPos < curMaxPos:
            self.adjCurPos[s] = curPos + 1
            if curPos + 1 < curMaxPos:
                self.nodesWUE[s] = len(self.path) - 1
            else:
                if s in self.nodesWUE:
                    del self.nodesWUE[s]
            v = self.adjList[s][curPos]
            self.path.append(v)
            s = v
            curPos = self.adjCurPos[s]
            curMaxPos = self.outDegree[s]
            self.numExploredEdges -= 1
        return

    def update_path(self, startPos):
        l = len(self.path) - 1
        self.path = self.path[startPos:l] + self.path[:startPos]
        for node, pos in self.nodesWUE.items():
            if pos < startPos:
                self.nodesWUE[node] = pos + l - startPos
            else:
                self.nodesWUE[node] = pos - startPos
        return

    def build_eulerian_cycle(self):
        self.explore(1)
        while self.numExploredEdges > 0:
            node, pos = self.nodesWUE.popitem()
            self.update_path(pos)
            self.explore(node)
        return self.path

    def print_path(self):
        print(' '.join([str(node + 1) for node in self.path[:-1]]))

    def read_input(self):
        data = list(sys.stdin.read().strip().split())
        self.numVertex, self.numExploredEdges = int(data[0]), int(data[1])
        self.adjList = [[] for _ in range(self.numVertex)]
        self.unusedEdges = [[] for _ in range(self.numVertex)]

        self.outDegree = [0] * self.numVertex
        self.inDegree = [0] * self.numVertex
        self.adjCurPos = [0] * self.numVertex

        for i in range(self.numExploredEdges):
            # find index of vertex
            curFrom = int(data[2 * i + 2]) - 1
            curTo = int(data[2 * i + 3]) - 1

            self.adjList[curFrom].append(curTo)
            self.outDegree[curFrom] += 1
            self.inDegree[curTo] += 1

        # Check In degree == Out degree
        for i in range(self.numVertex):
            if self.outDegree[i] != self.inDegree[i]:
                return False
        return True


if __name__ == "__main__":
    EulerianCycle()