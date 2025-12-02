#!/usr/bin/env python

# Originally written by James Torrance

import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
import os
import argparse

class Assembly:
     
    @property
    def original_zipped_fasta_file(self):
        return self._original_zipped_fasta_file

    @property
    def assembly_label(self):
        return self._assembly_label
    
    @property
    def fasta_suffix(self):
        return self._fasta_suffix
    
    @property
    def original_dir(self):
        return self._original_dir
    
    @property
    def zip_suffix(self):
        return self._zip_suffix
    
    def __init__(self, original_zipped_fasta_file):
        self._original_zipped_fasta_file = original_zipped_fasta_file
        (self._assembly_label, self._fasta_suffix, self._original_dir, self._zip_suffix) = self.get_parameters_from_filename(original_zipped_fasta_file)
        self.zipped_fasta_file = self.get_zipped_fasta_file(original_zipped_fasta_file)
	
    @staticmethod
    def get_parameters_from_filename(original_zipped_fasta_file):
        assembly_label = ''
        fasta_suffix = ''
        assembly_label_match = re.search('^([^/]+)\.(fa|fasta|fna)(\.gz)?$', os.path.basename(original_zipped_fasta_file))
        if assembly_label_match:
            original_dir = os.path.dirname(original_zipped_fasta_file)
            if original_dir != '':
                original_dir += '/'
            assembly_label = assembly_label_match.group(1)
            fasta_suffix = '.' + assembly_label_match.group(2)
            zip_suffix = assembly_label_match.group(3)
            if zip_suffix is None:
                zip_suffix = ''
        else:
            print('Error finding assembly label')
            quit()

        return (assembly_label, fasta_suffix, original_dir, zip_suffix)

    @staticmethod
    def get_zipped_fasta_file(original_zipped_fasta_file):
        zipped_fasta_file = ''
        fasta_file_match = re.search('^([^/]+)\.(fa|fasta|fna)(.gz)?$', os.path.basename(original_zipped_fasta_file))
        if fasta_file_match:
            zipped_fasta_file = fasta_file_match.group(1)
        else:
            print('Error finding fasta file')
            quit()

        return zipped_fasta_file

    def get_unzipped_screening_copy(self):
        unzipped_assembly_copy = self.original_dir + '/' + self.assembly_label + self.fasta_suffix

        # If the original was unzipped, just use that
        if self.zip_suffix == '':
            unzipped_assembly_copy = self.original_zipped_fasta_file

        return unzipped_assembly_copy

    def get_decontaminated_fasta_file(self):
        return(self.original_dir + self.assembly_label + '.decontaminated.fa')

    def get_contamination_file(self):
        contamination_file = self.original_dir + self.assembly_label + '.contamination'

        if os.path.isfile(contamination_file):
            return contamination_file
        else:
            return ''
        
    def get_bed_contamination_file(self):
        bed_contamination_file = os.path.join(self.original_dir, f"{self.assembly_label}.contamination.bed")
        return bed_contamination_file if os.path.isfile(bed_contamination_file) else ''
        
    def get_mito_fasta_file(self):
        return(self.original_dir + self.assembly_label + '.mito_removed.fa')

    def get_plastid_fasta_file(self):
        return(self.original_dir + self.assembly_label + '.plastid_removed.fa')

def main():

    parser = argparse.ArgumentParser(description='Remove contamination')
    parser.add_argument("--fasta", metavar='GRIT-42', type=str, help='FASTA files to clip')
    parser.add_argument("--bed", action='store_true', help='Use BED file')

    args = parser.parse_args()

    subassembly = Assembly(args.fasta)

    if args.bed:
        print(subassembly.get_unzipped_screening_copy(), subassembly.get_bed_contamination_file(), subassembly.get_decontaminated_fasta_file())
    else:
        print(subassembly.get_unzipped_screening_copy(), subassembly.get_contamination_file(), subassembly.get_decontaminated_fasta_file())
    
    verbose = True

    # INSTRUCTIONS: python3.4.0 clip_regions_oo.py [SUBASSEMBLY_LABEL]

    # NOTE: CHECK YOUR EXAMPLE DATA!

    coord_list_for_sequence = {}

    mito_flag = False
    plastid_flag = False
    common_euk_flag = False

    # MASK always converts to Ns
    # TRIM always trims (trim_Ns usually results in this)
    # CONTAMINANT converts to Ns unless it is at the very start or end, in which case it is trimmed (Vecscreen is CONTAMINANT)
    # MITOCHONDRIAL/PLASTID is like CONTAMINANT, but the sequence is set aside in a special mito file
    # REMOVE removes
    # Note that a REMOVE in the mitochondrial section results in both a REMOVE and a MITOCHONDRIAL entry
    # MITOCHONDRIAL_REMOVE or PLASTID_REMOVE directly signals the option given above

    contamination_file = subassembly.get_contamination_file()
    if args.bed:
        contamination_file = subassembly.get_bed_contamination_file()

    print(contamination_file)
    with open(contamination_file, 'r') as coord_input_handle:
        for line in coord_input_handle:
            if not re.match('^#', line):	# Ignore comments
                if args.bed:
                    fields = re.split('\s+', line)
                    id = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    if end == 0:
                        end = 1
                    treatment = fields[3]

                    if id not in coord_list_for_sequence:
                        coord_list_for_sequence[id] = []

                    coord_list_for_sequence[id].append([start+1,end,treatment])
                else:
                    if re.search('^\S*(REMOVE|TRIM|MASK|CLIP|CONTAMINANT|VecScreen)', line):
                        fields = re.split('\s+', line)
                        id = fields[1]

                        if re.search('^REMOVE', fields[0]):
                            treatment = 'REMOVE'
                            if id not in coord_list_for_sequence:
                                coord_list_for_sequence[id] = []
                            if mito_flag: # If this is mito sequence, flag it for copying too
                                coord_list_for_sequence[id].append([0,0, 'MITOCHONDRIAL_REMOVE'])
                            elif plastid_flag: # If this is plastid sequence, flag it for copying too
                                coord_list_for_sequence[id].append([0,0, 'PLASTID_REMOVE'])

                            else:
                                coord_list_for_sequence[id].append([0,0, treatment])

                        else:
                            start = fields[2]
                            end = fields[3]

                            treatment = 'MASK'
                            if re.search('^(TRIM|REVCLIP|FWDCLIP|CLIP)', fields[0]):
                                treatment = 'TRIM'
                            if re.search('^VecScreen', fields[0]):
                                treatment = 'CONTAMINANT'
                            if re.search('^CONTAMINANT', fields[0]):
                                treatment = 'CONTAMINANT'						

                            if id not in coord_list_for_sequence:
                                coord_list_for_sequence[id] = []
                            coord_list_for_sequence[id].append([int(start), int(end), treatment])
                    elif (mito_flag or common_euk_flag) and not re.search('^#', line) and not re.search('^\=', line):
                        fields = re.split('\s+', line)

                        if len(fields) > 6:
                            
                            treatment = 'CONTAMINANT'
                            if mito_flag:
                                treatment = 'MITOCHONDRIAL' # CAPTURE MITO
                            elif plastid_flag:
                                treatment = 'PLASTID' # CAPTURE MITO

                            id = fields[0]
                            start = int(fields[6])
                            end = int(fields[7])
                            start, end = sorted([start, end])
                            if id not in coord_list_for_sequence:
                                coord_list_for_sequence[id] = []
                            coord_list_for_sequence[id].append([start, end, treatment])
                    elif re.search('\=\=\=', line): # If this is a new section, this must be the end of the mito region
                        mito_flag = False
                        plastid_flag = False
                        common_euk_flag = False
                    if re.search('MITO', line): # If this is the mito region, change behaviour
                        mito_flag = True
                    elif re.search('PLASTID', line): # If this is the plastid region, change behaviour
                        plastid_flag = True
                    if re.search('COMMON CONTAMINANTS IN EUKARYOTES', line): # If this is the mito region, change behaviour
                        common_euk_flag = True

    # Merge entries
    #merged_coord_list_for_sequence = {}
    for id in coord_list_for_sequence:
        max_end = -100
        max_treatment = None
        for coord_pair_and_treatment in sorted(coord_list_for_sequence[id], key = lambda cpt: cpt[0]):
            if coord_pair_and_treatment[0] <= max_end:
                print('Screening coordinate overlap in ' + id + ' TREATMENT:' + coord_pair_and_treatment[2] + ' vs ' + max_treatment)
            if coord_pair_and_treatment[1] > max_end:
                max_end = coord_pair_and_treatment[1]
                max_treatment = coord_pair_and_treatment[2]

# NEED TO CHANGE MITO/REMOVE TO A SINGLE TYPE
# IF THERE IS A REGION OF OVERLAP, MERGE OR DIE
# IF MITO OVERLAPS- DIE
# IF SAME TYPE OVERLAPS- MERGE
# IF TRIM OVERLAPS CONTAMINANT- MERGE TO TRIM
# IF MASK OVERLAPS OTHER- DIE
# ANY OTHER COMBO- DIE
# PRINT RESULTS OF MERGING

        termini = []
        coord_pair_and_treatment_for_label = {}

        for coord_pair_and_treatment in sorted(coord_list_for_sequence[id], key = lambda cpt: cpt[0]):
        
            # Reject any malformed features
            if coord_pair_and_treatment[0] > coord_pair_and_treatment[1] or coord_pair_and_treatment[0] < 0 or coord_pair_and_treatment[1] < 0:
                print('Malformed feature:', str(coord_pair_and_treatment[0]), '-', str(coord_pair_and_treatment[1]))
                next

            label = coord_pair_and_treatment[2] + ':' + str(coord_pair_and_treatment[0]) + '-' + str(coord_pair_and_treatment[1])
            coord_pair_and_treatment_for_label[label] = coord_pair_and_treatment

            start_position = {
                'TERMINUS': 'START',
                'POSITION': coord_pair_and_treatment[0],
                'LABEL': label,
            }

            end_position = {
                'TERMINUS': 'END',
                'POSITION': coord_pair_and_treatment[1],
                'LABEL': label,
            }			

            termini.append(start_position)
            termini.append(end_position)
            
        termini = sorted(termini, key=sort_termini)

        depth = 0
        current_overlap_labels = []
        overlap_groups = []

        for terminus in termini:
            if terminus['TERMINUS'] == 'START':
                depth += 1
                current_overlap_labels.append(terminus['LABEL'])
            else:
                depth -= 1
                if depth == 0:
                    overlap_groups.append(current_overlap_labels)
                    current_overlap_labels = []


        merged_coord_pairs_and_treatments = []
        merge_flag = False

        for overlap_group in overlap_groups:
            #print('OVERLAP_GROUP')

            # If this label doesn't overlap, pass it to the merged set
            if len(overlap_group) == 1:
                merged_coord_pairs_and_treatments.append(coord_pair_and_treatment_for_label[overlap_group[0]])
            else:
                merge_flag = True
                starts = []
                ends = []
                treatments = {}
                for overlapping_label in overlap_group:
                    #print(overlapping_label, end=' ')
                    starts.append(coord_pair_and_treatment_for_label[overlapping_label][0])
                    ends.append(coord_pair_and_treatment_for_label[overlapping_label][1])
                    treatments[coord_pair_and_treatment_for_label[overlapping_label][2]] = 1

                starts = sorted(starts)
                ends = sorted(ends)
                min_start = starts[0]
                max_end = ends[-1]

                if len(treatments) == 1:
                    merged_coord_pair_and_treatment = [
                        min_start,
                        max_end,
                        list(treatments.keys())[0],
                    ]
                    merged_coord_pairs_and_treatments.append(merged_coord_pair_and_treatment)
                elif len(treatments) ==2 and 'TRIM' in treatments and 'CONTAMINANT' in treatments:
                    merged_coord_pair_and_treatment = [
                        min_start,
                        max_end,
                        'TRIM',
                    ]
                elif len(treatments) ==2 and 'REMOVE' in treatments and 'CONTAMINANT' in treatments:
                    merged_coord_pair_and_treatment = [
                        min_start,
                        max_end,
                        'REMOVE',
                    ]
                    merged_coord_pairs_and_treatments.append(merged_coord_pair_and_treatment)
                else:
                    quit('CANNOT MERGE: ' + str(overlap_group))

        if merge_flag:
            print('UNMERGED: ' + str(coord_list_for_sequence[id]))
            coord_list_for_sequence[id] = merged_coord_pairs_and_treatments
            print('  MERGED: ' + str(merged_coord_pairs_and_treatments))


    #coord_list_for_sequence = merged_coord_list_for_sequence

    mask_flag = False

    fasta_output_handle = open(subassembly.get_decontaminated_fasta_file(), 'w')

    mitochondrial_records = []
    plastid_records = []

    print('Opening', subassembly.get_unzipped_screening_copy())
    with open(subassembly.get_unzipped_screening_copy(), "r") as fasta_input_handle:
        for record in SeqIO.parse(fasta_input_handle, "fasta"):
            # print(record.id)

            remove_flag = False

            if(record.id in coord_list_for_sequence):
                if verbose:
                    print('Extracting from', record.id)

                last_end = 0
                edited_record = SeqRecord(Seq(''), id=record.id, description = record.description, name = record.name)
                for coord_pair_and_treatment in sorted(coord_list_for_sequence[record.id], key = lambda cpt: cpt[0]):
    
                    if verbose:
                        print('\t', coord_pair_and_treatment)

                    # We handle removal of mito sequence by treating it as a REMOVE entry and a MITOCHONDRIAL entry
                    # The MITOCHONDRIAL entry initially has coords -1,-1 until this stage, where we know the length
                    if (coord_pair_and_treatment[2] == 'MITOCHONDRIAL_REMOVE' or coord_pair_and_treatment[2] == 'PLASTID_REMOVE') and coord_pair_and_treatment[1] < 1:
                        coord_pair_and_treatment[0] = 1
                        coord_pair_and_treatment[1] = len(record.seq)

                    if coord_pair_and_treatment[2] == 'REMOVE' or coord_pair_and_treatment[2] == 'MITOCHONDRIAL_REMOVE' or coord_pair_and_treatment[2] == 'PLASTID_REMOVE':
                        remove_flag = True

                    if coord_pair_and_treatment[2] != 'REMOVE':
                        edited_record.seq += record.seq[last_end:(coord_pair_and_treatment[0]-1)] # Takes account of 0-based numbering
                        if coord_pair_and_treatment[2] == 'MASK' or ((coord_pair_and_treatment[2] == 'CONTAMINANT' or coord_pair_and_treatment[2] == 'MITOCHONDRIAL') and (coord_pair_and_treatment[0] > 1 and coord_pair_and_treatment[1] < len(record.seq))):
                            edited_record.seq += 'N' * (coord_pair_and_treatment[1] - coord_pair_and_treatment[0] + 1)
                        if coord_pair_and_treatment[2] == 'MITOCHONDRIAL' or coord_pair_and_treatment[2] == 'MITOCHONDRIAL_REMOVE':
                            mitochondrial_record = SeqRecord(Seq(''), id=record.id,description='', name = '')
                            mitochondrial_record.seq += record.seq[(coord_pair_and_treatment[0]-1):(coord_pair_and_treatment[1])]
                            mitochondrial_records.append(mitochondrial_record)
                        if coord_pair_and_treatment[2] == 'PLASTID' or coord_pair_and_treatment[2] == 'PLASTID_REMOVE':
                            plastid_record = SeqRecord(Seq(''), id=record.id,description='', name = '')
                            plastid_record.seq += record.seq[(coord_pair_and_treatment[0]-1):(coord_pair_and_treatment[1])]
                            plastid_records.append(plastid_record)

#							if coord_pair_and_treatment[1] < len(record.seq): # Create a new record for the remainder of the sequence if mito splits a record
#								SeqIO.write([edited_record], fasta_output_handle, 'fasta')
#								edited_record = SeqRecord(Seq(''), id=record.id + '_' + str(coord_pair_and_treatment[1]), description = record.description, name = record.name  + '_' + str(coord_pair_and_treatment[1]) )

                        last_end = coord_pair_and_treatment[1]

                edited_record.seq += record.seq[last_end:] # Finally add the end segment
                record = edited_record

            if not(remove_flag):
                record.id = re.sub('#', '_', record.id)
                record.description = re.sub('#', '_', record.description)
                record.name = re.sub('#', '_', record.name)

                SeqIO.write([record], fasta_output_handle, 'fasta')

    fasta_output_handle.close()

    if len(mitochondrial_records) > 0:
        mito_output_file = subassembly.get_mito_fasta_file()
        mito_output_handle = open(mito_output_file, 'w')
        SeqIO.write(mitochondrial_records, mito_output_handle, 'fasta')
        mito_output_handle.close()

    if len(plastid_records) > 0:
        plastid_output_file = subassembly.get_plastid_fasta_file()
        plastid_output_handle = open(plastid_output_file, 'w')
        SeqIO.write(plastid_records, plastid_output_handle, 'fasta')
        plastid_output_handle.close()

def sort_termini(terminus):
	if terminus['TERMINUS'] == 'END':
		terminus['POSITION'] += 0.5
	return terminus['POSITION']

if __name__ == "__main__":
	main()

# MERGE COORDS
# COLLECT ALL CONTAMINATION/FASTA FILES FROM DIRECTORIES (RE-USE OLD CODE)
