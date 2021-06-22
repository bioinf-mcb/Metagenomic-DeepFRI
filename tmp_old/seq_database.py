from bio_magic import sequence_similarity_score


import mmseqs





class SeqDatabase:
    def __init__(self):
        self.db_sequences = ['MIP', 'Swiss', 'PDB']
        self.db_distance_maps = ['MIP', 'Swiss', 'PDB']
        print("init DB")

    def find_similar(self, sequence, threshold):
        print("DB find_best_match")
        # in self.db_sequences
        output = []
        # for seq in self.db_sequences:
        if sequence_similarity_score("seq", "sequence") > threshold:
            output.append("seq")
        return output

    def get_distance_map(self, sequence):
        print("DB get distance map")
        # if sequence in self.db_sequences:
        #   return self.db_distance_maps[sequence]
        return None
