from Bio.Seq import Seq
#Creating Python Function To Convert DNA into Protein
def convert_dna_to_protein(dna_sequence):
    protein_sequence = Seq(dna_sequence).translate()
    return protein_sequence

#dna_sequence = "ATGCTGCCCGTAAAA"
#protein_sequence = convert_dna_to_protein(dna_sequence)
#print(f"Protein Sequence: {protein_sequence}")

#Crating Function to calculate hamming distance between slack username and twitter handle
def calculate_hamming_distance(slack_username, twitter_handle):
    hamming_distance = 0
    if len(slack_username) != len(twitter_handle):
        raise ValueError("Slack username and Twitter handle must have the same length")
    else:
        for i in range(len(slack_username)):
            if slack_username[i] != twitter_handle[i]:
                hamming_distance += 1
        return hamming_distance
#Example usage
slack_username = "@pavanpunnamraju"
twitter_handle = "pavan_punnamraju"
hamming_distance = calculate_hamming_distance(slack_username, twitter_handle)
print(f"Hamming Distance: {hamming_distance}")