"""Track how mutations spread geographically from large collections of
genomes.

WeightedTree:
    Assumptions:
        - No mutation reverses a previous mutation
        - Weights can be calculated relatively rapidly, but not so
            rapidly that every vertex could be checked against every
            vertex
    Algorithm:
        1. Start with a single vertex, then a single edge
        2. Add new vertex v by randomly choosing a leaf v_1 in the tree
            and finding the weight of a potential edge v-v_1
        3. Starting with v_i=v_1, Let v_{i+1} be the neighbor of v_i with
            the lowest weight of the potential edge v-v_{i+1} other than
            v_{i-1}
        4. When v-v_i has a lower weight than v-v_{i+1}, check to see if
            w(v,v_i) + w(v,v_{i+1} = w(v_i,v_{i+1} (v is intermediary to
            v_i and v_{i+1}). If so, subdivide the edge v_i-v_{i+1} to
            form the path v_i-v-v_{i+1}. If not, Add v as a leaf
            with edge v-v_i
        5. Check each vertex in N^2(v) (neighbors of neighbors of v) to
            see if they should be direct neighbors of v instead of
            second neighbors of v. Move them if so.

Weight Function:
    Assumptions:
        - Genomes are very similar
        - Any unknown reads "N" should be considered an exact match
            during comparisons
    Algorithm:
        1. Search for the next section that is identical.
        2. Extend the identical section as far as possible in both
            directions
        3. Apply Smith-Waterman between exact matches
        4. Return the sum of the number of mutations found by each
            application of Smith-Waterman

The process of calculating this tree requires a significant (hours) of
computer time. For this reason, the tree is saved as a file called
covid.tree after it is made. In successive uses of the program, if
covid.tree exists, then the tree is loaded instead of being recalculated
"""
import operator


class Genome:
    """Class to store a genome and accsession number"""

    def __init__(self, header_dict, genome_lines):
        self.accession = header_dict['Accession']
        self.header = header_dict
        genome_list = genome_lines.split("\n")
        self.genome = ''
        for line in genome_list:
            self.genome += line


    def get_genome(self, left="", right=""):
        """Get genome prepended and appended by left and right.

        Args:
            left (str): string to prepend with. (default = "")
            right (str): string to append with. (default = "")
        """
        return left + self.genome + right

    def make_fasta(self, header_list, width=70):
        """Make a fasta file from genome.

        Args:
            header_list (list): List of info to include in header
            width (int): Width to make each genome line (default=70)
        """
        out = ">"
        for key in header_list:
            if key in self.header:
                out += str(self.header[key]) + " "
        for ind in range(len(self.genome)):
            if ind % width == 0:
                out += "\n"
            out += self.genome[ind]
        return out


class WeightedTree:
    """A class to build trees one leaf at a time."""

    def __init__(self):
        self.verts = {}
        self.neighbors = {}
        self.edge_weights = {}
        self.root = None


    def get_edges(self):
        """ Returns a list of edges in the tree. """
        return list(self.edge_weights.keys())


    def get_leafs(self):
        """ Returns a list of leaves. """
        leaf_list = []
        for vert, neighbor_list in self.neighbors.items():
            if len(neighbor_list) == 1:
                leaf_list.append(vert)
        return leaf_list


    def reset_to_first_vertex(self, vertex_name, vertex_object):
        """Add the first vertex to the tree. This function is for private
        use only.
        """
        self.verts[vertex_name] = vertex_object
        self.neighbors = {vertex_name:[]}
        self.edge_weights = {}
        return True


    def save_tree(self, filename):
        """ Saves the tree to a file. Note: vertex objects are not saved"""
        with open(filename, 'w') as out_file:
            out_file.write("Neighbors\n")
            for vertex, neighbor_list in self.neighbors.items():
                out_file.write(vertex + ":")
                neighbors = ""
                for neighbor in neighbor_list:
                    neighbors += neighbor + ","
                out_file.write(neighbors[:-1] + "\n")
            out_file.write("Weights\n")
            for edge, weight in self.edge_weights.items():
                vert1, vert2 = edge
                out_file.write(vert1 + "," + vert2 + ":" + str(weight) + "\n")
            out_file.write("Root\n")
            if self.root is not None:
                out_file.write(self.root)
            else:
                out_file.write("None")


    def load_tree(self, filename, vertex_objects):
        """ Load a tree from file. Use objects from vertex_objects.

        Accession number is hard-coded
        """
        vert_list = []
        with open(filename) as in_file:
            line = in_file.readline()
            while line != "Weights\n":
                line = in_file.readline()
                if line != "Weights\n":
                    trimmed_line = line.rstrip()
                    vertex = trimmed_line.split(":")[0]
                    neighbors = trimmed_line.split(":")[1].split(",")
                    vert_list.append(vertex)
                    self.neighbors[vertex] = neighbors
            while line != "Root\n":
                line = in_file.readline()
                if line != "Root\n":
                    trimmed_line = line.rstrip()
                    weight = int(trimmed_line.split(":")[1])
                    edge = trimmed_line.split(":")[0].split(",")
                    self.edge_weights[(edge[0], edge[1])] = weight
            line = in_file.readline()
            if line != "None":
                self.root = line
        for vertex in vertex_objects:
            if vertex.accession in vert_list:
                self.verts[vertex.accession] = vertex


    def add_leaf(self, vertex_name, vertex_object, neighbor, weight):
        """Add a vertex to the tree.

        Args:
            name: unique identifier for the vertex (must be immutable)
            neighbor: name of neighbor
            weight: edge weight to neighbor
        """
        if len(self.verts) == 0:
            return self.reset_to_first_vertex(vertex_name, vertex_object)
        self.verts[vertex_name] = vertex_object
        self.neighbors[vertex_name] = [neighbor]
        self.neighbors[neighbor].append(vertex_name)
        self.edge_weights[edge_name(vertex_name, neighbor)] = weight
        return True


    def subdivide_edge(self, vertex_name, vertex_object, vert1, vert2, weight1,\
                       weight2):
        """Subdivide the edge v1 - v2 to form v1 - name - v2. edge
        weight of v1-name is weight1, edge weight of v2 - name is weight2
        """
        if weight1 + weight2 != self.edge_weights[edge_name(vert1, vert2)]:
            return False
        # Add vertex to vertex collection
        self.verts[vertex_name] = vertex_object
        # Remove vert1 -- vert2 edge
        self.neighbors[vert1].remove(vert2)
        self.neighbors[vert2].remove(vert1)
        self.edge_weights.pop(edge_name(vert1, vert2), None)
        # Add edge vert1 -- vertex_name edge
        self.neighbors[vert1].append(vertex_name)
        self.neighbors[vertex_name] = [vert1]
        self.edge_weights[edge_name(vert1, vertex_name)] = weight1
        # Add edge vert2 -- vertex_name edge
        self.neighbors[vert2].append(vertex_name)
        self.neighbors[vertex_name].append(vert2)
        self.edge_weights[edge_name(vert2, vertex_name)] = weight2
        return True


    def move_vertex(self, v_move_name, v_cur_name, v_new_name, new_weight):
        """ Move vertex v_move_name which is a neighbor of v_cur_name so
        that it becomes a neighbor of v_new_name with weight new_weight.
        """
        # Note: If v_move_name should have subdivided the edge v_cur_name - 
        # v_new_name, then original placement of v_move_name would have
        # been as an edge to v_new_name instead of to v_cur_name.
        self.neighbors[v_cur_name].remove(v_move_name)
        self.neighbors[v_move_name].remove(v_cur_name)
        self.neighbors[v_new_name].append(v_move_name)
        self.neighbors[v_move_name].append(v_new_name)
        self.edge_weights.pop(edge_name(v_cur_name, v_move_name))
        self.edge_weights[edge_name(v_move_name, v_new_name)] = new_weight
        return True


    def smallest_neighbor(self, tree_vertex_name, \
                          tree_vertex_wt, new_vertex_obj, weight_func):
        """Returns (neighbor, distance) of tree_vertex that is closest to
        new_vertex as measured by weight_func. Assumes that at most one
        neighbor will be closer than tree_vertex.)

        Args:
            used_vertex_name: vertex identifier of vertex that has
                already been used (Keep vertices from bouncing back and
                forth)
            tree_vertex_name: vertex identifier of vertex whose neigbors
                will be checked
            tree_vertex_wt (numeric): Weight of potential edge from
                new_vertex to tree_vertex
            new_vertex_obj: object of the vertex to be added
            weight_func (function): A function object used to evaluate
                the weight of potential edges.
        """
        smallest_vert = None
        smallest_wt = None
        for vertex in self.neighbors[tree_vertex_name]:
            cur_wt = weight_func(self.verts[vertex], new_vertex_obj)
            if cur_wt is None:
                return None
            if cur_wt < tree_vertex_wt:
                return (vertex, cur_wt)  # At most 1 vertex smaller
            if smallest_wt is None or cur_wt < smallest_wt:
                smallest_wt = cur_wt
                smallest_vert = vertex
        return smallest_vert, smallest_wt


    def check_neighbor_move(self, new_vertex, old_vertex, weight_func):
        """ Check all the neighbors of old_vertex (except new_vertex) to
        see if they should be moved to being neighbors or new_vertex. If
        so, then move them.
        """
        for vert in self.neighbors[old_vertex]:
            if vert != new_vertex:
                new_weight = weight_func(self.verts[vert],\
                                         self.verts[new_vertex])
                if new_weight < self.edge_weights[edge_name(vert, old_vertex)]:
                    self.move_vertex(vert, old_vertex, new_vertex, new_weight)


    def add_vertex(self, vertex_name, vertex_obj, weight_func):
        """Add a vertex (leaf or subdivide an edge) to its nearest
        neighbor based on the calculated value in weight_func

        Args:
            vertex_name: The identifier for a vertex
            vertex_obj: The object that the vertex represents
            weight_func (function): A function object used to evaluate
                the weight of potential edges.
        """
        if vertex_name in self.verts:
            return False
        if len(self.verts) == 0:
            return self.reset_to_first_vertex(vertex_name, vertex_obj)
        if len(self.verts) == 1:
            cur_vertex_name = list(self.verts.keys())[0]
            cur_weight = weight_func(vertex_obj, self.verts[cur_vertex_name])
            if cur_weight is None:
                return None
            return self.add_leaf(vertex_name, vertex_obj,\
                                 cur_vertex_name, cur_weight)
        decreasing = True
        cur_vertex_name, cur_vertex_obj = list(self.verts.items())[0]
        cur_weight = weight_func(cur_vertex_obj, vertex_obj)
        if cur_weight is None:
            return None
        while decreasing:
            try:
                next_vertex_name, next_vertex_wt = \
                    self.smallest_neighbor(cur_vertex_name, cur_weight, \
                                           vertex_obj, weight_func)
            except TypeError:
                return None
            if self.subdivide_edge(vertex_name, vertex_obj, cur_vertex_name,\
                        next_vertex_name, cur_weight, next_vertex_wt):
                self.check_neighbor_move(vertex_name, cur_vertex_name, weight_func)
                self.check_neighbor_move(vertex_name, next_vertex_name, weight_func)
                return True
            if next_vertex_wt >= cur_weight:
                decreasing = False
            else:
                cur_weight = next_vertex_wt
                cur_vertex_name = next_vertex_name
                cur_vertex_obj = self.verts[cur_vertex_name]
        return_val = self.add_leaf(vertex_name, vertex_obj, cur_vertex_name,\
                cur_weight)
        self.check_neighbor_move(vertex_name, cur_vertex_name, weight_func)
        return return_val


    def set_root(self, vertex_name):
        """Set the root vertex"""
        self.root = vertex_name


    def output_tree_str(self):
        """Return a string containing a breadth first or depth first
        textual representation of the graph.
        """
        if self.root is None:
            start = list(self.verts.keys())[0]
        else:
            start = self.root
        return self.depth_first_out_str(start, None)


    def depth_first_out_str(self, start, pred, sep=","):
        """Return a depth first representation of the tree rooted at
        start.
        """
        out_str = ''
        out_deg = 0
        for neighbor in self.neighbors[start]:
            if neighbor != pred:
                out_deg += 1
        if out_deg == 0:
            out_str += str(start) + "\n"
            return out_str
        for neighbor in self.neighbors[start]:
            if neighbor != pred:
                out_str += str(start) + sep + \
                        self.depth_first_out_str(neighbor, start)
        return out_str


    def vert_to_vert(self, vertex1, vertex2):
        """Returns a path from vertex1 to vertex2."""
        if vertex1 == vertex2:
            return [vertex1]
        distance = {0:[vertex2]}
        predecessor = {0:[]}
        vertex1_found = False
        next_dist = 1
        while not vertex1_found:
            distance[next_dist] = []
            predecessor[next_dist] = []
            for vert in distance[next_dist - 1]:
                for nbr in self.neighbors[vert]:
                    if nbr not in predecessor[next_dist - 1]:
                        distance[next_dist].append(nbr)
                        predecessor[next_dist].append(vert)
                        if nbr == vertex1:
                            vertex1_found = True
            next_dist += 1
        next_dist -= 1  # Go back to farthest distance
        path = []
        cur_vert = vertex1
        while next_dist > 0:
            path.append(cur_vert)
            cur_vert = \
                    predecessor[next_dist][distance[next_dist].index(cur_vert)]
            next_dist -= 1
        path.append(cur_vert)
        return path


    def path_to_root(self, leaf):
        """Returns a list of vertices from a leaf back to a vertex."""
        if self.root is None:
            root = list(self.verts.keys())[0]
        else:
            root = self.root
        return self.vert_to_vert(leaf, root)


    def path_to_leaf(self, leaf):
        """Returns a list of vertices from the root to a leaf."""
        path = self.path_to_root(leaf)
        if isinstance(path, list):
            path.reverse()
        return path


    def output_adjacency_entry(self, vert_i, vert_j):
        """ Return the i,j-entry of the adjacency matrix """
        if edge_name(vert_i, vert_j) in self.edge_weights.keys():
            return 1
        return 0


    def save_adjacency_as_csv(self, filename, include_names=True):
        """ Write a csv file for the adjacency matrix.

        This can be used for programs like graphviz to visualize the
        tree.
        """
        with open(filename, 'w') as outfile:
            line = ""
            if include_names:
                for vert in self.verts:
                    line += "," + str(vert)
                outfile.write(line + "\n")
                line = ""
            for verti in self.verts:
                if include_names:
                    line = str(verti) + ","
                else:
                    line = ""
                for vertj in self.verts:
                    line += str(self.output_adjacency_entry(verti, vertj)) + ","
                outfile.write(line[:-1] + "\n")


def get_accession_from_fasta_header(header):
    """ Take a fasta header line and return the accession number."""
    header = header[1:]  # Remove leading >
    header = header.split("|")[0]  # First portion
    header = header.rstrip()
    if '.' in header:  # Remove version
        header = header.split('.')[0]
    return header


def fasta_to_genomes(filename, read_info, filters=[('Length', '>', 29000)]):
    """ Returns a list of Genome objects from the fasta file given by filename.
    Filters based on the tuple filter_tuple using data from read_info.
    Args:
        filename (str): Name of fasta file.
        read_info (dict of dicts): Dictionary of genome information
            collected from csv file
        filter_tuple (tuple): tuple of key, conditional, condition where
            key is the key from read_info[current] to filter on
            conditional is one of <, <=, !=, ==, >=, >
            condition is the condition to check.
    """
    genome_list = []
    cur_header = ''
    cur_fasta = ''
    with open(filename) as fasta:
        for line in fasta:
            if line[0] == ">":  # header line
                if cur_fasta != "":  # fasta data
                    accession = get_accession_from_fasta_header(cur_header)
                    if filter_by_dict(read_info[accession], filters):
                        genome_list.append(Genome(read_info[accession], \
                                                  cur_fasta))
                    cur_fasta = ""
                cur_header = line.rstrip()
            else:
                cur_fasta += line
        genome_list.append(Genome(read_info[accession], cur_fasta))
    return genome_list


def filter_by_dict(data_info, filters):
    """ Returns True if the dictionary data_info meets the condition in
    filter_tuple.

    Args:
        data_info (dict): Dictionary of overview data
        filter_tuple (tuple): tuple of key, conditional, condition where
            key is the key from read_info[current] to filter on
            conditional is one of <, <=, ==, >=, >, in
            condition is the condition to check.
    """
    for filter_tuple in filters:
        key, conditional, condition = filter_tuple
        if conditional == '<':
            oper = operator.lt
        elif conditional == '<=':
            oper = operator.le
        elif conditional == '==':
            oper = operator.eq
        elif conditional == '!=':
            oper = operator.ne
        elif conditional == '>=':
            oper = operator.ge
        elif conditional == '>':
            oper = operator.gt
        elif conditional == 'in':
            oper = operator.contains
        try:
            if not oper(data_info[key], condition):
                return False
        except TypeError:
            if '.' in data_info[key]:
                if not oper(float(data_info[key]), condition):
                    return False
            else:
                if not oper(int(data_info[key]), condition):
                    return False
    return True


def edge_name(vert1, vert2):
    """Sort vert1 and vert2 (names) and return as a tuple to be used
    for edge names
    """
    if vert1 < vert2:
        return (vert1, vert2)
    return (vert2, vert1)


def csv_line_to_dict(head_list, line, sep=','):
    """ Return a dictionary from a line from a CSV file. Assumes that
    anything in quotes should not be separated.

    Args:
        head_list (list of strings): List of headings to be used as keys
        line (str): Line from the CSV file
        sep (str): separator (, by default)
    """
    line = line.rstrip()
    out_dict = {}
    line = line.replace(',"', '"')
    line = line.replace('",', '"')
    line = line.split('"')
    entries = []
    for index, entry in enumerate(line):
        if index % 2 == 0:
            entries.extend(entry.split(sep))
        else:
            entries.append(entry)
    for index, entry in enumerate(entries):
        out_dict[head_list[index]] = entry
    return out_dict


def read_csv_file_to_dict(filename, key_col, sep=','):
    """ Take a CSV file and return a dictionary where the keys are given
    by the key_col entry. Assumes that the first line of the file is a
    header that defines the meaning of each column.

    Args:
        filename (str): Name of csv file to read
        key_col (str): The column to use as keys in the dictionary.
        sep (str): The character separating columns
    """
    csv_dict = {}
    try:
        with open(filename) as csv_file:
            header_line = next(csv_file)
            header_line = header_line.rstrip()  # Remove trailing \n
            header = header_line.split(sep)
            counter = 0
            for line in csv_file:
                line_dict = csv_line_to_dict(header, line.rstrip(), sep)
                key = line_dict.get(key_col, str(counter))
                counter += 1
                csv_dict[key] = line_dict
        return csv_dict
    except FileNotFoundError:
        return {}


def matches(strand1, strand2):
    """ Return True if strand1 and strand2 are identical."""
    pairs = zip(strand1, strand2)
    for (char1, char2) in pairs:
        if char1 == char2:
            pass
        elif char1 == 'N' or char2 == 'N':
            pass
        else:
            return False
    return True


def find_next_identity(start, strand1, strand2, min_match=50, max_offset=50):
    """Find the next area where strand1 and strand2 exactly match in at
    least min_match spaces. Expand that until the match is maximal.

    Args:
        start (int, int): Indices for (strand1, strand2) where the search
            will start from.
        strand1 (str): A strand of DNA
        strand2 (str): A strand of DNA
        min_match (int): The minimum size of a match to consider.
        max_offset (int): Maximum distance to search down strand1.
    """
    start1, start2 = start
    try:
        min_strand2 = \
                strand2[start2 + max_offset: start2 + max_offset + min_match]
        check_strand1 = \
                strand1[start1: start1 + 3 * max_offset]
    except IndexError:
        return None
    count = 1
    while min_strand2 not in check_strand1:
        count += 1  # Grow the window to check
        try:
            min_strand2 = \
                strand2[start2 + count * max_offset:\
                        start2 + count * max_offset + min_match]
            #check_strand1 = \
            #    strand1[start1 + (count - 1) * max_offset:\
            #            start1 + (count + 2) * max_offset]
            check_strand1 = \
                strand1[start1 + max_offset:\
                        start1 + (count + 2) * max_offset]
        except IndexError:
            return None
    index1 = strand1.find(min_strand2, start1 + (count - 1) * max_offset, \
             start1 + (count + 2) * max_offset)
    index2 = start2 + count * max_offset
    length = min_match
    # Decrease
    while matches(strand1[index1:index1 + length], \
                  strand2[index2:index2 + length]) and \
                        index1 >= 0 and index2 >= 0:
        index1 -= 1
        index2 -= 1
        length += 1
    # Went 1 too far
    index1 += 1
    index2 += 1
    length -= 1
    # Increase
    try:
        while matches(strand1[index1 + length], \
                strand2[index2 + length]):
            length += 1
        # Went 1 too far
        length -= 1
    except IndexError: #  Reached end of strand
        length -= 1
    return (index1, index2, length)


def compute_matrices(seq1, seq2, cost_dict):
    """Compute a mutation cost and a mutation path from genome1 to
    genome2

    Args:
        seq1 (str): First genome
        seq2 (str): Second genome
        cost_dict (dict): dictionary of costs for mutations
    """
    # make predecessor and cost matrices
    pred_matrix = []
    cost_matrix = []

    # Make first row of pred_matrix and cost_matrix
    pred_row = ['left'] * len(seq1)
    pred_row[0] = 'start'
    cost_row = [0]
    for ind in range(len(seq1) - 1):
        cost_row.append(cost_row[ind] + cost_dict['insert'])
    pred_matrix.append(pred_row)
    cost_matrix.append(cost_row)

    # Make remaining rows
    for ind2 in range(1, len(seq2)):
        pred_row = ['up']
        cost_row = [cost_dict['delete']]
        for ind1 in range(1, len(seq1)):
            left = cost_dict['insert']
            above = cost_dict['delete']
            diag = cost_dict['change']
            if seq1[ind1] == 'N' or seq2[ind2] == 'N' or\
                    seq2[ind2] == seq1[ind1]:
                diag = 0
            left = cost_row[ind1 - 1] + left
            above = cost_matrix[ind2 - 1][ind1] + above
            diag = cost_matrix[ind2 - 1][ind1 - 1] + diag
            min_cost = min(left, above, diag)
            cost_row.append(min_cost)
            if diag == min_cost:
                pred_row.append('diag')
            elif above == min_cost:
                pred_row.append('up')
            else:
                pred_row.append('left')
        pred_matrix.append(pred_row)
        cost_matrix.append(cost_row)
    return cost_matrix, pred_matrix


def compute_cost(genome1, genome2, cost_dict=\
        {'insert':1, 'delete':1, 'change':1}, max_matrix=5000\
        ):
    """Read a predecessor matrix and return the path from start to
    end.
    Args:
        genome2 (Genome): Second genome
        genome1 (Genome): First genome
        cost_dict (dict): dictionary of costs for mutations
    """
    seq1 = genome1.get_genome()
    seq2 = genome2.get_genome()
    index1 = 0
    index2 = 0
    cost = 0
    done = False
    while not done:
        try:
            (new_index1, new_index2, match_length) = \
                    find_next_identity((index1, index2), seq1, seq2)
        except TypeError:  # No matches found til end of sequence
            skip_section1 = "N" + seq1[index1:] + "N"
            skip_section2 = "N" + seq2[index2:] + "N"
            matrices = compute_matrices(skip_section1,\
                                        skip_section2, cost_dict)
            cost += matrices[0][-1][-1]
            done = True
        if new_index1 - index1 >= max_matrix or\
                new_index2 - index2 >= max_matrix:
            return None
        if new_index1 == index1 and new_index2 == index2:  # no skips
            index1 += match_length
            index2 += match_length
        elif new_index1 == index1:  # skip in strand2 only
            cost += new_index2 - index2
            index2 = new_index2 + match_length
            index1 += match_length
        elif new_index2 == index2:  # skip in strand1 only
            cost += new_index1 - index1
            index1 = new_index1 + match_length
            index2 += match_length
        else:  # skip in both (change?)
            skip_section1 = "N" + seq1[index1:new_index1] + "N"
            skip_section2 = "N" + seq2[index2:new_index2] + "N"
            matrices = compute_matrices(skip_section1,\
                                        skip_section2, cost_dict)
            cost += matrices[0][-1][-1]
            index1 = new_index1 + match_length
            index2 = new_index2 + match_length
        if index1 >= len(seq1) - 2 and index2 >= len(seq2) - 2:  # done
            done = True
        elif index1 >= len(seq1) - 2:  # done with seq2
            cost += len(seq2) - 2 - index2
            done = True
        elif index2 >= len(seq2) - 2:  # done with seq1
            cost += len(seq1) - 2 - index1
            done = True
    return cost


def compute_path(genome1, genome2, cost_dict, show_pred=False):
    """Read a predecessor matrix and return the path from start to
    end.

    Args:
        genome1 (Genome): First genome
        genome2 (Genome): Second genome
        cost_dict (dict): dictionary of costs for mutations
    """
    cost_matrix, pred_matrix = compute_matrices(genome1, genome2, cost_dict)
    if show_pred:
        seq1 = genome1.get_genome("N", "N")
        seq2 = genome2.get_genome("N", "N")
        print("  " + seq1)
        for ind in range(len(pred_matrix)):
            print(seq2[ind], end=" ")
            for base in pred_matrix[ind]:
                if base == 'left':
                    print("-", end="")
                elif base == 'up':
                    print("|", end="")
                elif base == 'diag':
                    print("\\", end="")
                elif base == 'start':
                    print("*", end="")
            print()
    y_index = len(pred_matrix) - 1
    x_index = len(pred_matrix[y_index]) - 1
    cur_read = pred_matrix[y_index][x_index]
    backwards_path = []
    direction_reverse = {'up': 'down', 'left': 'right', 'diag': 'diag'}
    # Read path backwards
    while cur_read != 'start':
        backwards_path.append(direction_reverse[cur_read])
        if cur_read == 'diag':
            x_index -= 1
            y_index -= 1
        elif cur_read == 'left':
            x_index -= 1
        else:
            y_index -= 1
        cur_read = pred_matrix[y_index][x_index]

    # Reverse path
    path = []
    while len(backwards_path) > 0:
        path.append(backwards_path.pop())
    return cost_matrix[-1][-1], path


def initialize_state_dict():
    """ Creates a dictionary using all 50 states where
    dict['AL']['NY'] represents a weight of an edge from Alabama to NY.
    Set all weights to 0 to be updated later.
    """
    states = ["AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL",\
            "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD",\
            "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ",\
            "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",\
            "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY", "PR", "DC"]
    state_dict = {}
    for state1 in states:
        state_dict[state1] = {state2:0 for state2 in states}
    return state_dict


def get_state(location):
    """ Return a state abbreviation from location info. """
    states = {"AK":"Alaska", "AL":"Alabama", "AR":"Arkansas", "AZ":"Arizona",\
            "CA":"California", "CO":"Colorado", "CT":"Connecticut",\
            "DE":"Delaware", "FL":"Florida", "GA":"Georgia", "HI":"Hawaii",\
            "IA":"Iowa", "ID":"Idaho", "IL":"Illinois", "IN":"Indiana",\
            "KS":"Kansas", "KY":"Kentucky", "LA":"Louisiana",\
            "MA":"Massachusetts", "MD":"Maryland", "ME":"Maine",\
            "MI":"Michigan", "MN":"Minnesota", "MO":"Missouri",\
            "MS":"Mississippi", "MT":"Montana", "NC":"North Carolina",\
            "ND":"North Dakota", "NE":"Nebraska", "NH":"New Hampshire",\
            "NJ":"New Jersey", "NM":"New Mexico", "NV":"Nevada",\
            "NY":"New York", "OH":"Ohio", "OK":"Oklahoma", "OR":"Oregon",\
            "PA":"Pennsylvania", "RI":"Rhode Island", "SC":"South Carolina",\
            "SD":"South Dakota", "TN":"Tennessee", "TX":"Texas", "UT":"Utah",\
            "VA":"Virgina", "VT":"Vermont", "WA":"Washington",\
            "WI":"Wisconsin", "WV":"West Virginia", "WY":"Wyoming",\
            "PR":"Puerto Rico", "DC":"Washington DC"}
    for code, state in states.items():
        if code.upper() in location.upper():
            return code
        if state.upper() in location.upper():
            return code
    return None


def earlier(str_date1, str_date2):
    """ Returns str_date1 if str_date1 is earlies than str_date2 and vice
    versa. Returns None if they are the same day.
    """
    if str_date1 == str_date2:
        return None
    try:
        year1 = int(str_date1.split("-")[0])
        month1 = int(str_date1.split("-")[1])
        day1 = int(str_date1.split("-")[2])
        year2 = int(str_date2.split("-")[0])
        month2 = int(str_date2.split("-")[1])
        day2 = int(str_date2.split("-")[2])
    except IndexError:
        return None
    # Check year
    if year1 < year2:
        return str_date1
    if year2 < year1:
        return str_date2
    # Check month
    if month1 < month2:
        return str_date1
    if month2 < month1:
        return str_date2
    # Check day
    if day1 < day2:
        return str_date1
    if day2 < day1:
        return str_date2
    return None


def find_oldest(metadata, accession_list):
    """ Return a list of accession number of all the sequences that have
    the oldest date using data from metadata.
    """
    oldest_date = metadata[accession_list[0]]["Collection_Date"]
    oldest_accession = []  # accession added in for loop
    for accession in accession_list:
        date = metadata[accession]["Collection_Date"]
        older = earlier(oldest_date, date)
        if older is None:
            oldest_accession.append(accession)
        elif older == date:
            oldest_date = date
            oldest_accession = [accession]
    return oldest_accession


def find_extreme_dates(metadata, accession_list):
    """ Returns (oldest, newest) date in accession_list
    """
    oldest_date = metadata[accession_list[0]]["Collection_Date"]
    newest_date = metadata[accession_list[0]]["Collection_Date"]
    if len(oldest_date.split('-')) == 2:
        oldest_date += "-01"
    if len(newest_date.split('-')) == 2:
        newest_date += "-31"
    for accession in accession_list:
        date = metadata[accession]["Collection_Date"]
        older = earlier(oldest_date, date)
        newer = earlier(newest_date, date)
        if older == date:
            oldest_date = date
        if newer == newest_date:
            newest_date = date
    return (oldest_date, newest_date)


def fan_out_from_accession_list(accession_list, tree, metadata):
    """ Return list of edges in Tree that fan out from the vertices in
    accession_list.
    """
    leaf_list = tree.get_leafs()
    edge_list = []
    length = len(leaf_list)
    num = 0
    for leaf in leaf_list:
        num += 1
        print("\r" + "Finding path from " + str(num) + " of " +\
                str(length) + " leafs. " + str(len(accession_list)) + " verts in accession list.", end="")
        path_to_root = tree.path_to_root(leaf)
        # Determine if the path to the root goes through a vertex in
        # the accession list so will be used.
        in_accession_list = False
        for accession in accession_list:
            if accession in path_to_root:
                in_accession_list = True
                break
        if in_accession_list:
            # Add to the accession list until you reach the first
            # vertex that is already in the accession list
            new_accession = True
            for vert1, vert2 in zip(path_to_root, path_to_root[1:]):
                if new_accession:
                    accession_list.append(vert1)
                    edge_list.append(edge_name(vert1, vert2))
                    if vert2 in accession_list:
                        #print(" ",metadata[vert1]["Geo_Location"],metadata[vert2]["Geo_Location"],edge_name(vert1,vert2) in edge_list, end=" ")
                        new_accession = False
                        break
    print()
    return edge_list


def next_date(date):
    """ Returns the day after date as a string in YYYY-MM-DD format. """
    def feb_days(year):
        """ Returns the number of days in Feb in that year. """
        if year % 4 == 0:
            return 29
        return 28
    def int_to_two(num):
        """ Returns a two digit string of the integer (<100) num with
        a leading zero if necessary.
        """
        if num < 10:
            return "0" + str(num)
        return str(num)
    try:
        year = int(date.split("-")[0])
        month = int(date.split("-")[1])
        day = int(date.split("-")[2])
    except KeyError:
        return None
    days_mo = {1:31, 2:feb_days(year), 3:31, 4:30, 5:31, 6:30, 7:31,
               8:31, 9:30, 10:31, 11:30, 12:31}
    day += 1
    if day > days_mo[month]:
        day = 1
        month += 1
    if month > 12:
        month = 1
        year += 1
    return str(year) + "-" + int_to_two(month) + "-" + int_to_two(day)


def days_list(date1, date2):
    """ Returns a list of days between date1 and date2 inclusive.
    Thus, if date1 == date2, then returns 1
    Dates are assumed to be in YYYY-MM-DD format.
    """
    older = earlier(date1, date2)
    if older is None:  # Same day doesn't need to process further
        return [date1]
    old_date = older.strip()
    new_date = date2.strip()
    if old_date == new_date:
        new_date = date1
    cur_date = old_date
    dates = [old_date]
    while cur_date != new_date:
        if cur_date is not None:
            cur_date = next_date(cur_date)
            dates.append(cur_date)
        else:
            cur_date = new_date
    return dates


def edges_output(filename, edges, metadata):
    """ Writes a CSV file showing number of edges from each state to
    each state in the list edges.
    """
    states_dict = initialize_state_dict()
    for edge in edges:
        accession1, accession2 = edge
        state1 = get_state(metadata[accession1]["Geo_Location"])
        state2 = get_state(metadata[accession2]["Geo_Location"])
        date1 = metadata[accession1]["Collection_Date"]
        date2 = metadata[accession2]["Collection_Date"]
        if state1 is None or state2 is None:
            if state1 is None:
                print(metadata[accession1]["Geo_Location"])
            if state2 is None:
                print(metadata[accession2]["Geo_Location"])
        else:
            older = earlier(date1, date2)
            if older is None:
                if state1 == state2:
                    states_dict[state1][state2] += 1
                else:
                    pass
            if date1 == older:
                states_dict[state1][state2] += 1
            else:
                states_dict[state2][state1] += 1
    states = ["AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL",\
            "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD",\
            "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ",\
            "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",\
            "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY"]
    with open(filename, 'w') as out_csv:
        out_line = ','  # skip first column
        for state in states:
            out_line += state + ','
        out_csv.write(out_line[:-1])  # skip trailing comma
        out_csv.write("\n")
        for state in states:
            out_line = state + ","
            for state2 in states:
                out_line += str(states_dict[state][state2]) + ","
            out_csv.write(out_line[:-1])
            out_csv.write("\n")


def make_animation(accession_list, tree, metadata, base_filename):
    """ Make a series of csv files holding data from
    fan_out_from_accession_list.
    """
    print("Root:", metadata[tree.root]["Geo_Location"],
          metadata[tree.root]["Collection_Date"])
    edges = fan_out_from_accession_list(accession_list, tree, metadata)
    all_verts = []
    for vert1, vert2 in edges:
        if vert1 not in all_verts:
            all_verts.append(vert1)
        if vert2 not in all_verts:
            all_verts.append(vert2)
    geo_list = []
    for vert in all_verts:
        if get_state(metadata[vert]["Geo_Location"]) not in geo_list:
            geo_list.append(get_state(metadata[vert]["Geo_Location"]))
    for geo in geo_list:
        print(geo,end=", ")
    start_date, end_date = find_extreme_dates(metadata, all_verts)
    print("Start Date:", start_date, "End Date:", end_date)
    dates = {}
    cur_date = start_date
    while cur_date != end_date:
        dates[cur_date] = initialize_state_dict()
        cur_date = next_date(cur_date)
    dates[end_date] = initialize_state_dict()
    print("Edges found: " + str(len(edges)))
    skipped = 0
    MA_count = 0
    for edge in edges:
        accession1, accession2 = edge
        state1 = get_state(metadata[accession1]["Geo_Location"])
        state2 = get_state(metadata[accession2]["Geo_Location"])
        if state2 == "USA: MA" or state2 == "USA: MA":
            MA_count += 1
        date1 = metadata[accession1]["Collection_Date"]
        date2 = metadata[accession2]["Collection_Date"]
        if state1 is None or state2 is None:
            if state1 is None:
                print("State Error:", accession1,\
                      metadata[accession1]["Geo_Location"])
                skipped += 1
            if state2 is None:
                print("State Error:", accession2,\
                      metadata[accession2]["Geo_Location"])
                skipped += 1
        elif date1 is None or date2 is None:
            if date1 is None:
                print("Date Error:", accession1,\
                      metadata[accession1]["Collection_Date"])
                skipped += 1
            if date2 is None:
                print("Date Error:", accession1,\
                      metadata[accession2]["Collection_Date"])
                skipped += 1
        else:
            #print(state1 + ", " + date1 + " - " + state2 + ", " + date2)
            date_list = days_list(date1, date2)
            num_days = len(date_list)
            for day in date_list:
                try:
                    dates[day][state1][state2] += 1 / num_days
                except KeyError:
                    skipped += 1
    print(str(MA_count) + " edges have endpoint in MA.", end="")
    print(str(skipped) + " skipped")
    states = ["AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL",\
            "GA", "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD",\
            "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ",\
            "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN",\
            "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY"]
    for date in dates.keys():
        with open(base_filename + "-" + date + ".csv", 'w') as out_csv:
            out_line = ','  # skip first column
            for state in states:
                out_line += state + ','
            out_csv.write(out_line[:-1])  # skip trailing comma
            out_csv.write("\n")
            for state in states:
                out_line = state + ","
                for state2 in states:
                    out_line += str(dates[date][state][state2]) + ","
                out_csv.write(out_line[:-1])
                out_csv.write("\n")


def main():
    filters_list = [('Length', '>', 29000),
                    ('Nuc_Completeness', '==', 'complete'),
                    ('Geo_Location', 'in', 'USA:')]
    metadata = read_csv_file_to_dict('sequences.csv', 'Accession')
    genomes = fasta_to_genomes('sequences.fasta', metadata, filters_list)
    skipped = []
    print(str(len(genomes)) + " selected out of " + \
          str(len(list(metadata.keys()))) + " genomes in file.")

    filename = "covid.tree"
    load = True
    try:
        temp = open(filename)
        temp.close()
    except FileNotFoundError:
        load = False

    genome_tree = WeightedTree()
    if load:
        genome_tree.load_tree(filename, genomes)
    else:
        count = 0
        for genome in genomes:
            count += 1
            if genome_tree.add_vertex(genome.accession, genome,
                                      compute_cost) is None:
                skipped.append(genome.accession)
            else:
                print("\r" + str(count + 1) + ' done (' + \
                      str(len(skipped)) + " skipped)")
        genome_tree.save_tree(filename)
    edges = genome_tree.get_edges()
    print("Number of vertices in tree: " + str(len(edges)+1))
    edges_output("states_spread.csv", edges, metadata)
    genome_tree.set_root(find_oldest(metadata,\
                         list(genome_tree.verts.keys()))[0])
    massachusetts_all = []
    for vert in genome_tree.verts.keys():
        if get_state(metadata[vert]["Geo_Location"]) == "MA":
            massachusetts_all.append(vert)
    #mass_start = find_oldest(metadata, massachusetts_all)
    make_animation(massachusetts_all, genome_tree, metadata, "MA_spread/MA_spread")


if __name__ == "__main__":
    main()
