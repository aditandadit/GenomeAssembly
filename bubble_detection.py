# python3
import itertools
import sys

'''
1) The De Bruijn graph has hanging tips, which means that some nodes don't have any outgoing edge.
2) There can be more than one bubble between two nodes. The PDF isn't completely clear about this.
3) The size of a bubble is determined by the number of edges, not number of vertices.
'''

'''
https://en.wikipedia.org/wiki/DNA_read_errors#Tips_and_Bubbles
Given nodes 𝑣 and 𝑤, a (𝑣, 𝑤)-bubble => a pair of short non-overlapping disjoint paths between 𝑣 and 𝑤.
Bubble is formed when an error occurs during sequence reading and 
there is a path for k-mer reads to reconnect with main graph

Try and look at what bubbles look like => ABCD is a bubble
    A----B
---       ---
    C----D

'''


class bubble_detection:

    def __init__(self, k_mer_size, threshold, reads):
        self.k = k_mer_size  # break Reads into k-mers of size k
        self.threshold = threshold  # Threshold above which bubble length needs to be to be detected

        self.de_bruijn_graph = {}  # De bruijn graph
        self.paths = {}  # path during dgs

        self.num_bubbles = 0
        # for each node, we're storing [set of outgoing edges, number of incoming edges]
        self.build_de_bruijn_graph(self.reads_to_kmers(reads))

    def reads_to_kmers(self, reads):
        k_mers = []
        for read in reads:
            for j in range(len(read) - self.k + 1):
                k_mers.append(read[j: j + self.k])
        return k_mers

    def build_de_bruijn_graph(self, kmers):
        """
            How to Build de-bruijn graph:
            1. For each k_mer get left and right K-1 mer
            2. Add an edge from left -> right in the graph
            3. If left and right already exist in the graph use those nodes as opposed to creating new ones
        """
        for kmer in kmers:
            left, right = kmer[:-1], kmer[1:]

            if left != right:
                # for each node, add [set of outgoing edges, number of incoming edges] by default
                self.de_bruijn_graph.setdefault(left, [set(), 0])
                self.de_bruijn_graph.setdefault(right, [set(), 0])

                if right not in self.de_bruijn_graph[left][0]:
                    self.de_bruijn_graph[left][0].add(right)
                    self.de_bruijn_graph[right][1] += 1

    def count_bubbles(self):
        for k, v in self.de_bruijn_graph.items():
            # If node has multiple outgoing edges its can be start of bubble
            if self.has_multiple_outgoing_edges(k):
                self.dfs(path=[k], start=k, current=k, depth=0)

        # paths has candidates of [... [V --> W]...]
        # where V has multiple outgoing, and W has multiple incoming and bubble length is beyond a threshold
        # for Combination of these paths where their intersection size = 2 => Bubble
        for _, candidates_list in self.paths.items():
            # get every pair in paths -> nC2
            for pair in itertools.combinations(candidates_list, r=2):
                # # only v,w ie the start and end of bubble are shared
                if len(set(pair[0]) & set(pair[1])) == 2:
                    self.num_bubbles += 1

        return self.num_bubbles

    def dfs(self, path, start, current, depth):
        # start -> V, current -> W
        # Bubble -> Non overlapping paths between V and W
        # start always has multiple outgoing, current always has multiple incoming
        if current != start and self.has_multiple_incoming(current):
            self.paths \
                .setdefault((start, current), list()) \
                .append(path[:])

        if depth == self.threshold:
            return

        # continue dfs for neighbors of current
        for neighbor in self.de_bruijn_graph[current][0]:
            if neighbor not in path:  # avoid cycle
                path.append(neighbor)
                self.dfs(path, start, neighbor, depth + 1)
                path.remove(neighbor)

    def has_multiple_incoming(self, vertex):
        return self.de_bruijn_graph[vertex][1] > 1

    def has_multiple_outgoing_edges(self, vertex):
        return len(self.de_bruijn_graph[vertex][0]) > 1


if __name__ == "__main__":
    data = sys.stdin.read().split()
    k, t, reads = data[0], data[1], data[2:]
    print(bubble_detection(int(k), int(t), reads).count_bubbles())
