# python3
'''
Task.
Break the reads into fragments of 𝑘 = 15
Before constructing the de Bruijn graph, remove tips and handle bubbles.

Dataset.
Each line of the input contains a single read.
Each read is 100 nucleotides long and contains a single sequencing error
in order to simulate the 1% error rate of Illumina sequencing machines.
You are not given the 100-mer composition of the genome
'''

import sys
from collections import deque


class DeBruijnGraph(object):

    def __init__(self, k, reads):
        self.k = k
        self.threshold = self.k + 1
        self.kmer_id_2Way_map = Kmer_Id_2Way_Map()
        self.coverage = {}
        self.de_bruijn_graph = {}

        self.num_outgoing = lambda k: len(self.de_bruijn_graph[k][0])
        self.num_incoming = lambda k: self.de_bruijn_graph[k][1]

        self.build_de_bruijn_graph(self.reads_to_kmers(reads))

    def reads_to_kmers(self, reads):
        break_read = lambda read: [read[j:j + self.k] for j in range(len(read) - self.k + 1)]
        return [kmer for read in reads for kmer in break_read(read)]

    def build_de_bruijn_graph(self, kmers):
        for kmer in kmers:
            left = self.kmer_id_2Way_map.insert(kmer[:-1])
            right = self.kmer_id_2Way_map.insert(kmer[1:])

            if left != right:
                self.de_bruijn_graph.setdefault(left, [set(), 0])
                self.de_bruijn_graph.setdefault(right, [set(), 0])
                self.coverage.setdefault((left, right), 0)
                self.coverage[(left, right)] += 1

                if right not in self.de_bruijn_graph[left][0]:
                    self.de_bruijn_graph[left][0].add(right)
                    self.de_bruijn_graph[right][1] += 1

    def remove_leaves(self):
        removable = [k for k, v in self.de_bruijn_graph.items() if len(v[0]) == 0]
        for k in removable:
            del self.de_bruijn_graph[k]

    def print_graph(self):
        for k, v in self.de_bruijn_graph.items():
            print(k, v)


class RemoveTips(DeBruijnGraph):

    def __init__(self, k, reads):
        DeBruijnGraph.__init__(self, k, reads)

    def remove_tips(self):
        for k, v in self.de_bruijn_graph.items():
            find_and_remove = None

            if self.num_outgoing(k) == 1 and self.num_incoming(k) == 0:
                find_and_remove = self.remove_inward
            elif self.num_outgoing(k) > 1:
                find_and_remove = self.remove_outward
            else:
                continue

            condition = True
            while condition:
                condition = False
                for edge in v[0]:
                    if find_and_remove(edge, 0):
                        v[0].remove(edge)
                        condition = True
                        break

    def remove_outward(self, current, depth):
        if self.num_outgoing(current) > 1 or self.num_incoming(current) > 1:
            return False

        if depth == self.threshold:
            return False

        if self.num_outgoing(current) == 0:
            return True

        if self.remove_outward(next(iter(self.de_bruijn_graph[current][0])), depth + 1):
            to = next(iter(self.de_bruijn_graph[current][0]))
            self.de_bruijn_graph[current][0].pop()
            self.de_bruijn_graph[to][1] -= 1
            return True

        return False

    def remove_inward(self, current, depth):
        if self.num_outgoing(current) == 0 or self.num_incoming(current) > 1:
            return True

        if depth == self.threshold:
            return False

        if self.remove_inward(next(iter(self.de_bruijn_graph[current][0])), depth + 1):
            to = next(iter(self.de_bruijn_graph[current][0]))
            self.de_bruijn_graph[current][0].pop()
            self.de_bruijn_graph[to][1] -= 1
            return True

        return False


class RemoveBubbles(RemoveTips):

    def __init__(self, k, reads):
        RemoveTips.__init__(self, k, reads)
        self.paths = {}

    def remove_bubbles(self):
        for k, v in self.de_bruijn_graph.items():
            if self.num_outgoing(k) > 1:
                self.dfs(path=[k], current=k, depth=0)

        for pair, candidates_list in self.paths.items():
            source, target = pair[0], pair[1]
            best_path = max(candidates_list, key=lambda item: item[1])[0]
            for path, _ in candidates_list:
                if best_path == path or not self.bubble_possible(source, target):
                    continue
                if self.paths_disjoint(best_path, path) and self.path_exists(path):
                    self.remove_path(path)

    def bubble_possible(self, source, target):
        return len(self.de_bruijn_graph[source][0]) > 1 and self.de_bruijn_graph[target][1] > 1

    def remove_path(self, path):
        for j in range(len(path) - 1):
            self.de_bruijn_graph[path[j]][0].remove(path[j + 1])
            self.de_bruijn_graph[path[j + 1]][1] -= 1
            del self.coverage[(path[j], path[j + 1])]

    def paths_disjoint(self, a, b):
        return len(set(a) & set(b)) == 2

    def path_exists(self, path):
        for j in range(len(path) - 1):
            if path[j + 1] not in self.de_bruijn_graph[path[j]][0]:
                return False
        return True

    def dfs(self, path, current, depth):
        if current != path[0] and self.num_incoming(current) > 1:
            weight = sum(self.coverage[(path[i], path[i + 1])] for i in range(len(path) - 1)) / len(path)
            self.paths.setdefault((path[0], current), list()).append((path[:], weight))

        if depth == self.threshold:
            return

        for next_ in self.de_bruijn_graph[current][0]:
            if next_ not in path:
                path.append(next_)
                self.dfs(path, next_, depth + 1)
                path.remove(next_)


class PhiX174GenomeAssembler(RemoveBubbles):

    def __init__(self, k, reads):
        RemoveBubbles.__init__(self, k, reads)

    def make_eulerian_cycle(self):
        vertices = deque()
        path = []
        current = next(iter(self.de_bruijn_graph))
        vertices.append(current)

        while vertices:
            current = vertices[0]
            if len(self.de_bruijn_graph[current][0]) != 0:
                t = next(iter(self.de_bruijn_graph[current][0]))
                vertices.append(t)
                self.de_bruijn_graph[current][0].remove(t)
                continue
            path.append(current)
            vertices.popleft()

        return path

    def assemble_genome(self):
        self.remove_tips()
        self.remove_leaves()
        self.remove_bubbles()

        cycle = self.make_eulerian_cycle()
        circular_genome = self.kmer_id_2Way_map.id_to_kmer_map[cycle[0]]
        for i in range(1, len(cycle) - (self.k - 1)):
            circular_genome += self.kmer_id_2Way_map.id_to_kmer_map[cycle[i]][-1]

        return circular_genome


class Kmer_Id_2Way_Map:
    def __init__(self):
        self.num_ids = 0
        self.kmer_to_id_map = {}
        self.id_to_kmer_map = {}

    def insert(self, kmer):
        if kmer not in self.kmer_to_id_map:
            self.kmer_to_id_map[kmer] = self.num_ids
            self.id_to_kmer_map[self.num_ids] = kmer
            self.num_ids += 1
        return self.kmer_to_id_map[kmer]


if __name__ == "__main__":
    print(PhiX174GenomeAssembler(20, sys.stdin.read().split()).assemble_genome())
