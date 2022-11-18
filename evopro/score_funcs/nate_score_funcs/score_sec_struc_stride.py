def stride_formatted(pdb):
    import subprocess
    import os
    
   # try:
   #    subprocess.check_output(["stride", f"{pdb}", "-o"], stderr=subprocess.STDOUT)
   # except subprocess.CalledProcessError as e:
   #     print(e.output)
   # else:
   #     string_out = subprocess.check_output(["stride", f"{pdb}", "-o"], stderr=subprocess.STDOUT)
   #     string_out_split = string_out.decode().split("\n")
    string_out = subprocess.check_output(["stride", f"{pdb}", "-o"], stderr=subprocess.STDOUT)
    string_out_split = string_out.decode().split("\n")

#runs stride and breaks up output

    all_structure = []

    for line in string_out_split:
        if line.startswith('LOC'):
            col = line.split(" ")
            l = [x for x in col if x]
            all_structure.append(l)

    helix_strand_structure  = []
    for i in all_structure:
        if i[1] == 'AlphaHelix':
            x = ['H', int(i[3]), int(i[6])]
            helix_strand_structure.append(x)
        if i[1] == 'Strand':
            y = ['E', int(i[3]), int(i[6])]
            helix_strand_structure.append(y)
    helix_strand_structure.sort(key=lambda x: x[1])
#creates list of [Structure, Start, End] sorted by start position
    
    corrected_structure = []
    
    for i in range(len(helix_strand_structure)-1):
        if helix_strand_structure[i][0] == helix_strand_structure[i+1][0]:
            if int(helix_strand_structure[i+1][1])-int(helix_strand_structure[i][2]) <= 1:
                corrected_structure.append([helix_strand_structure[i][0], helix_strand_structure[i][1], 
                    helix_strand_structure[i+1][2]])
                i+=1
            else:
                corrected_structure.append(helix_strand_structure[i])
        else:
            corrected_structure.append(helix_strand_structure[i])
#combines [Structure, Start, End] for two structures with one res gap. Appends original structure to list otherwise.
    if corrected_structure:
        if helix_strand_structure[-1][2] != corrected_structure[-1][2]:
            corrected_structure.append(helix_strand_structure[-1])
#makes sure last structure is included in list

    return corrected_structure
   
def align_structures(target, pdb):
    from needlemanwunsch import nw
    import numpy as np
    
    target_list = []
    af2_list = []

    for i in target:
        target_list.append(i[0])
    for i in pdb:
        af2_list.append(i[0])
    
    target_string = ''.join(target_list)
    af2_string = ''.join(af2_list)
    alignment = nw(target_string, af2_string)
    alignment_split = alignment.split("\n")
    target_string_aligned = alignment_split[0]
    af2_string_aligned = alignment_split[1]

    return target_string_aligned, af2_string_aligned
#creates string for target and af2 with aligned secondary structure blocks

def insert_gaps(stride_output, string_aligned):
    j = 0
    updated_output = []
    gap = ['-', 0, 0]

    for t, i in zip(string_aligned, range(len(string_aligned))):
        if t  != '-':
            updated_output.append(stride_output[j])
            j += 1
        else:
            updated_output.append(gap)

    return updated_output
#inserts gaps into corrected, aligned stride outputs

def determine_score(final_target_output, final_af2_output):
    structure_score = 0
    length_score = 0

    for i in range(len(final_target_output)):
        if final_target_output[i][0] == final_af2_output[i][0]:
            structure_score += 1
            target_structure_length = final_target_output[i][2]-final_target_output[i][1]+1
            af2_structure_length = final_af2_output[i][2]-final_af2_output[i][1]+1
            length_score +=  (abs(target_structure_length-af2_structure_length)/target_structure_length)
        if final_target_output[i][0] != final_af2_output[i][0]:
            length_score += 1

    percent_structure_score = (structure_score/len(final_target_output))*100
    average_length_score = (length_score/len(final_target_output))*100   

    return percent_structure_score, average_length_score
#average_structure_score is percent of secondary structure blocks that are correct
#average_length_score is average percent error in length for each correct block
#Note that this function assumes the target has at least one secondary structure block we want
#to recreate.  If none of the af2 blocks match the target blocks, the length score is zero.
#Both scores are percentages (two digit numbers, not decimals)

def score_sec_struc(path_to_target, pdb):
    import subprocess
    from needlemanwunsch import nw
    import numpy as np
    import os   

    with open('temp.pdb', 'w') as f:
        f.write(pdb)

    target_formatted = stride_formatted(path_to_target)
    af2_formatted = stride_formatted('temp.pdb')
    target_aligned, af2_aligned = align_structures(target_formatted, af2_formatted)
    target_gaps = insert_gaps(target_formatted, target_aligned)
    af2_gaps = insert_gaps(af2_formatted, af2_aligned)
    structure_score, length_score = determine_score(target_gaps, af2_gaps)

    return structure_score, length_score

#path_to_target = 'ferr_ems_00120.pdb'
#path_to_pdb = 'HHH_bc_16535.pdb'
#a, b = score_sec_struc(path_to_target, path_to_pdb)
#print(a)
#print(b)

