# Program checks to see if sequences share overlapping bases between their 5' and 3' ends.


class FastaObj:
    def __init__(self, fasta_id, sequence):
        self.fasta_id = fasta_id
        self.sequence = sequence

    def get_id(self):
        return self.fasta_id

    def get_seq(self):
        return self.sequence


def create_fasta_list(file_name):
    # Get FASTA ID and corresponding sequence for each sequence in text file.
    file = open(file_name, "r")

    current_id = ""
    current_seq = ""
    fasta_list = []

    for line in file:
        # New FASTAs are marked with a ">". Everytime we encounter a ">" we want to take the ID which is the remaining
        # characters on the first line, then concatenate all subsequent lines until the next ">" as the sequence
        # for that ID. Whenever a new FASTA is discovered the previous one is stored in a list, the final FASTA is
        # is added on completion of the loop.
        if line[0] == ">" and current_id == "":
            current_id = line[1:].rstrip()
        elif line[0] == ">" and current_id != "":
            fasta_list.append(FastaObj(current_id, current_seq))
            current_id = line[1:].rstrip()
            current_seq = ""
        else:
            current_seq += line.rstrip()
    else:
        fasta_list.append(FastaObj(current_id, current_seq))

    file.close()
    return fasta_list


def find_overlap(input_list, overlap_num):
    # Parameters of function are a list of FASTAs and the number of bases to check for overlap.
    overlap_list = []
    for fasta_item in input_list:
        item_seq = fasta_item.get_seq()
        item_id = fasta_item.get_id()
        for check_item in input_list:
            # Check each FASTA against all other FASTAs.
            check_seq = check_item.get_seq()
            check_id = check_item.get_id()
            if item_seq[-overlap_num:] == check_seq[:overlap_num] and item_id != check_id:
                # Check if the 5' end of the first FASTA matches the 3' end of the second FASTA.
                directed_edge = [item_id, check_id]
                # print(directed_edge)
                overlap_list.append(directed_edge)

    return overlap_list


object_list = create_fasta_list("rosalind_grph.txt")
edge_list = find_overlap(object_list, 3)
print(edge_list)

output_file = open("overlap_output.txt", "w")
for line in edge_list:
    output_line = line[0] + " " + line[1] + "\n"
    output_file.write(output_line)
output_file.close()
