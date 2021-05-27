from seq_database import SeqDatabase
from docker_communication import get_input_sequence, output_functionality
from deep_fri_core import DeepLM, DeepFRI, DeepCNN
from bio_magic import pairwise_alignment, distance_map_quality


# constant
seq_similarity_threshold = 0.95
contact_threshold = 6
dm_quality_threshold = 0.9

new_sequence = get_input_sequence()

predicted_functionality = None

# try to find distance map necessary for DeepFRI_GCN
distance_map = None
# find multiple similar sequences across databases
db = SeqDatabase()
similar_sequences = db.find_similar(new_sequence, seq_similarity_threshold)

if len(similar_sequences) > 0:
    # check quality of multiple distance maps across similar sequences
    best_dm_score = -1
    for seq in similar_sequences:
        dm = db.get_distance_map(seq)
        dm = pairwise_alignment(dm, seq, new_sequence)
        dm_score = distance_map_quality(dm)
        if dm_score > dm_quality_threshold and dm_score > best_dm_score:
            best_dm_score = dm_score
            distance_map = dm

# generate distance map on the fly using language model
if distance_map is None:
    deep_lm = DeepLM()
    lm_distance_map = deep_lm(new_sequence)
    if distance_map_quality(lm_distance_map) > dm_quality_threshold:
        distance_map = lm_distance_map

# run GCN using obtained distance map
if distance_map is not None:
    deep_fri = DeepFRI()
    contact_map = distance_map < contact_threshold
    predicted_functionality = deep_fri(new_sequence, contact_map)

# use CNN if satisfying distance map was not obtained
else:
    deep_cnn = DeepCNN()
    predicted_functionality = deep_cnn(new_sequence)

output_functionality(predicted_functionality)
