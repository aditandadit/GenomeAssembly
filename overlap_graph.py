# Uses python3

DEFAULT_READS_NUMBER = 1618  # num reads
DEFAULT_MIN_OVERLAP_LENGTH = 70  # minimum overlap needed
LENGTH_OF_READ = 100  # number of bases in read


class TrieNode(object):
    def __init__(self):
        # children map
        self.children = {}
        # index from reads array that this node represents
        self.indexes = []


# Trie of Prefix of Strings, But only above certain length
class PrefixTrie(object):
    def __init__(self):
        self.root = TrieNode()

    def add(self, read, index):
        # For Prefixes of length min_overlap and above only
        for end in range(DEFAULT_MIN_OVERLAP_LENGTH, len(read)):
            rev_prefix = read[:end][::-1]

            # Add prefix into trie
            self.add_to_trie(index, rev_prefix)

    def add_to_trie(self, index, rev_prefix):
        node = self.root
        for char in rev_prefix:
            if char not in node.children:
                node.children[char] = TrieNode()
            node = node.children[char]

        # add index of the read to node
        node.indexes.append(index)

    def match(self, string):
        """
        Traverse the Trie and match strings above default overlap threshold
        :return list of match tuples (index, length)
        """
        adjacent = []
        node = self.root
        length = 0

        for char in string[::-1]:
            if char not in node.children:
                break
            node = node.children[char]
            length += 1

            if length >= DEFAULT_MIN_OVERLAP_LENGTH and node.indexes:
                for index in node.indexes:
                    adjacent.append((index, length))
        return adjacent


def overlap_value_bw_strings(s, t):
    for i in range(LENGTH_OF_READ, 0, -1):
        if s[LENGTH_OF_READ - i:] == t[:i]:
            return i
    return 0


def build_overlap_graph(reads):
    prefixTrie = PrefixTrie()
    for index, read in enumerate(reads):
        prefixTrie.add(read, index)

    overlap_graph = [[] for _ in range(len(reads))]
    for index, read in enumerate(reads):
        # get Reads that are neighbors to this read from prefix trie
        overlap_graph[index] = prefixTrie.match(read)
        # sort the neighbors based on overlap length
        overlap_graph[index].sort(key=lambda neighbor: neighbor[1], reverse=True)

    return overlap_graph


def build_hamiltonian_path_greedy(adj):
    current = 0
    added = {0}
    path = [(0, 0)]
    while len(added) < len(adj):
        for i, link in enumerate(adj[current]):
            if link[0] not in added:
                added.add(link[0])
                current = link[0]
                path.append(link)
                break
    return path


def assemble_genome(hamiltonian_path, reads):
    genome = ""
    for node in hamiltonian_path:
        genome += reads[node[0]][node[1]:]
    genome = genome[:-overlap_value_bw_strings(reads[hamiltonian_path[-1][0]], reads[0])]
    return genome


def read_inputs():
    reads = []
    for i in range(DEFAULT_READS_NUMBER):
        reads.append(input())
    reads = list(set(reads))
    return reads

'''
Step 1: Build Overlap Graph
Step 2: Build Hamiltonian Path in Greedy Fashion 
(This doesn't get optimal solution as Hamiltonian Path has no polynomial time solution)
'''
reads = read_inputs()
overlap_graph = build_overlap_graph(reads)
hamiltonian_path = build_hamiltonian_path_greedy(overlap_graph)

genome = assemble_genome(hamiltonian_path, reads)
print(genome)
