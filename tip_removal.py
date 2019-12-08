# python3
import sys

'''
Tips are error-prone ends of the reads that do not form a bubble but instead form a path starting in a vertex
without incoming edges or ending in a vertex without outgoing edges in the de Bruijn graph.
Tips should be removed iteratively because removing a tip can expose another tip
'''

'''
Task.
* Given a list of error-prone reads,
* Construct a de Bruijn graph from the 15-mers created from the reads
* Perform the task of tip removal on this de Bruijn graph.

Dataset:
400 reads of length 100, each on a separate line.
Each read is 100 nucleotides long and contains a single sequencing error (i.e., one mismatch per read)
in order to simulate the 1% error rate of Illumina sequencing machines.
Note that you are not given the 100-mer composition of the genome
(i.e., some 100-mers may be missing).
'''


class remove_tip:

    def __init__(self, k, reads):
        self.k = k
        self.threshold = self.k
        self.de_bruijn_graph = {}
        self.paths = {}
        self.edges_removed = 0
        # for each node, we're storing [set of outgoing edges, number of incoming edges]
        self.build_de_bruijn_graph(self.reads_to_kmers(reads))
        # print(self.de_bruijn_graph)

    def reads_to_kmers(self, reads):
        k_mers = []
        for read in reads:
            for j in range(len(read) - self.k + 1):
                k_mers.append(read[j: j + self.k])
        # print(k_mers)
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


    def remove_tips(self):
        """
        Iterate over nodes, for nodes which could be tips run a dfs and keep removing tips
        Removing a Tip also exposes Other tips so take care of that
        """
        for node, edge_info in self.de_bruijn_graph.items():
            find_and_remove = None
            # 2 cases of remove
            # case 1 => Single Outgoing edge, no incoming edge
            # case 2 => more than 1 outgoing edge
            if len(edge_info[0]) == 1 and edge_info[1] == 0:
                find_and_remove = self.remove_incoming
            elif len(edge_info[0]) > 1:
                find_and_remove = self.remove_outgoing
            else:
                continue

            condition = True
            while condition:
                condition = False
                for edge in edge_info[0]:
                    if find_and_remove(edge, 0):
                        edge_info[0].remove(edge)
                        self.edges_removed += 1
                        condition = True
                        break

        return self.edges_removed

    def remove_outgoing(self, current, depth):
        if self.num_outgoing(current) > 1 or self.num_incoming(current) > 1:
            return False  # This is definitely not a tip

        if self.num_outgoing(current) == 0:
            return True  # this is a tip

        if depth == self.threshold:
            return False

        # Recurse and if it returns true => We removed a tip and now another tip is exposed
        if self.remove_outgoing(next(iter(self.de_bruijn_graph[current][0])), depth + 1):
            self.de_bruijn_graph[current][0].pop()
            self.edges_removed += 1
            return True

        return False

    def remove_incoming(self, current, depth):
        if depth == self.threshold:
            return False

        if self.num_outgoing(current) == 0 or self.num_incoming(current) > 1:
            return True

        if self.remove_incoming(next(iter(self.de_bruijn_graph[current][0])), depth + 1):
            self.de_bruijn_graph[current][0].pop()
            self.edges_removed += 1
            return True

        return False

    def num_incoming(self, v):
        return self.de_bruijn_graph[v][1]

    def num_outgoing(self, v):
        return len(self.de_bruijn_graph[v][0])


if __name__ == "__main__":
    k_mer_size = 15
    k, reads = k_mer_size, sys.stdin.read().split()

    print(remove_tip(k, reads).remove_tips())

'''
Another Interesting Approach (Think about it later)
Find Strongly Connected Components, Tip will never be part of any SCC
Keep Doing till untill all tips are gone?
'''
