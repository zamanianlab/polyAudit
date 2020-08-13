import sys
import argparse
import pathlib
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import csv
import pandas as pd


def primertrim(records, primers, output):
    # remove PacBio primers
    trimmed_records = []
    with open(primers) as primer_file:
        primer_list = []
        for line in primer_file:
            primer_list.append(line.strip())

    with open(records) as fasta:
        records_iter = SeqIO.parse(fasta, 'fasta')
        for record in records_iter:
            for primer in primer_list:
                len_primer = len(primer)
                if record.seq.endswith(primer):
                    # slice off primer
                    record.seq = record.seq[:-len_primer]
                    trimmed_records.append(record)
                    break
                else:
                    pass
    output_path = pathlib.Path(output)
    input_path = pathlib.Path(records)
    stem = input_path.stem
    with open(pathlib.Path.joinpath(output_path, stem + "_trimmed.fasta"), 'w') as trimmed:
        for record in trimmed_records:
            SeqIO.write(record, trimmed, 'fasta')
    return trimmed_records


def measurepolya(trimmed_records, output_directory, input_file):
    '''
    Measure the polyA tails of trimmed sequences.
    '''
    ids = []
    tails = []
    if type(trimmed_records) is list:
        print("Using trimmed sequences to measure poly(A) tails.")
        for record in trimmed_records:
            index = 1
            t_counter = 1
            mismatch_counter = 0
            rc = record.seq.reverse_complement()
            for nt in rc:
                if nt == "T":
                    t_counter += 1
                    index += 1
                    mismatch_counter = 0
                else:
                    index += 1
                    mismatch_counter += 1
                    if mismatch_counter == 2:
                        ids.append(record.id)
                        tails.append(t_counter)
                        break
                    else:
                        pass
    elif type(trimmed_records) is str:
        print("Using provided sequences to measure poly(A) tails.")
        with open(trimmed_records) as fasta:
            records_iter = SeqIO.parse(fasta, 'fasta')
            for record in records_iter:
                index = 1
                t_counter = 1
                mismatch_counter = 0
                rc = record.seq.reverse_complement()
                for nt in rc:
                    if nt == "T":
                        t_counter += 1
                        index += 1
                        mismatch_counter = 0
                    else:
                        index += 1
                        mismatch_counter += 1
                        if mismatch_counter == 2:
                            ids.append(record.id)
                            tails.append(t_counter)
                            break
                        else:
                            pass
    dict = {ids[i]: tails[i] for i in range(len(ids))}
    output_path = pathlib.Path(output_directory)
    input_path = pathlib.Path(input_file)
    stem = input_path.stem
    with open(pathlib.Path.joinpath(output_path, stem + "_polya.csv"), 'w') as polyas:
        for key in dict.keys():
            polyas.write("%s,%s\n" % (key, dict[key]))
    return dict


def findkmer(trimmed_records, length_file, output_directory, k, n, input_file):
    '''
    Removes poly(A) tails and then searches the upstream sequence for
    overrepresented kmers. kmer length provided by -k, which defaults to 6.
    Returns a df with the top 10 (non AAAAAA or TTTTTT) k-mers.
    '''
    counts = {}
    if type(length_file) is dict:
        lengths = length_file
    if type(length_file) is str:
        with open(length_file) as input:
            reader = csv.reader(input)
            lengths = {rows[0]: rows[1] for rows in reader}
    if type(trimmed_records) is list:
        for record in trimmed_records:
            try:
                length = int(lengths.get(record.id))
            except TypeError:
                print("Length cannot be coerced to type int.")
                continue
            polya = int(lengths.get(record.id))
            start = len(record.seq) - (50 + polya)
            sliced_seq = str(record.seq[start:-polya])
            num_kmers = len(sliced_seq) - k + 1
            for i in range(num_kmers):
                # Slice the string to get the kmer
                kmer = sliced_seq[i:i + k]
                # Add the kmer to the dictionary if it's not there
                if kmer not in counts:
                    counts[kmer] = 0
                    # Increment the count for this kmer
                    counts[kmer] += 1
                else:
                    counts[kmer] += 1
    elif type(trimmed_records) is str:
        with open(trimmed_records) as fasta:
            records_iter = SeqIO.parse(fasta, 'fasta')
            for record in records_iter:
                try:
                    length = int(lengths.get(record.id))
                except TypeError:
                    print("Length cannot be coerced to type int.")
                    continue
                polya = int(lengths.get(record.id))
                start = len(record.seq) - (50 + polya)
                sliced_seq = str(record.seq[start:-polya])
                num_kmers = len(sliced_seq) - k + 1
                for i in range(num_kmers):
                    # Slice the string to get the kmer
                    kmer = sliced_seq[i:i + k]
                    # Add the kmer to the dictionary if it's not there
                    if kmer not in counts:
                        counts[kmer] = 0
                        # Increment the count for this kmer
                        counts[kmer] += 1
                    else:
                        counts[kmer] += 1

    out = pd.DataFrame.from_dict(counts, orient='index',
                                 columns=['Count'])
    out.index.name = 'kmer'
    out.reset_index(inplace=True)
    out.sort_values(by=['Count'], inplace=True, ascending=False)
    output_path = pathlib.Path(output_directory)
    input_path = pathlib.Path(input_file)
    stem = input_path.stem
    with open(pathlib.Path.joinpath(output_path, stem + "_" + str(k) + "mer.csv"), 'w') as kmers:
        out.to_csv(kmers, index=False)
    out = out[out.kmer != 'AAAAAA']
    out = out[out.kmer != 'TTTTTT']
    out = out[:n]
    return out


def findpas(trimmed_records, length_file, kmer_file, output_directory, input_file):
    '''
    Removes poly(A) tails and searches the upstream sequences for kmers
    provided by the user or piped from kmercount(). Starts with the top kmer,
    so will only find non-canonical PAS motifs if there hasn't been a canonical
    (AAUAAA) already found.
    '''
    records = []
    seqs = []
    indices = []
    found = []

    kmers = []
    if isinstance(kmer_file, pd.DataFrame):
        kmers = kmer_file['kmer'].tolist()
    elif type(kmer_file) is str:
        with open(kmer_file) as input:
            for line in input:
                kmers.append(line.strip())
    # check which type of lengths were provided
    if type(length_file) is dict:
        lengths = length_file
    if type(length_file) is str:
        with open(length_file) as input:
            reader = csv.reader(input)
            lengths = {rows[0]: rows[1] for rows in reader}
    if type(trimmed_records) is list:
        for record in trimmed_records:
            polya = int(lengths.get(record.id))
            start = len(record.seq) - (50 + polya)
            sliced_seq = str(record.seq[start:-polya])
            for kmer in kmers:
                if kmer in sliced_seq:
                    index = sliced_seq.find(kmer)
                    index = index - 5
                    records.append(record.id)
                    out_seq = str(record.seq[start:-(polya - 15)])
                    seqs.append(out_seq)
                    indices.append(index)
                    found.append(kmer)
                    break
                elif kmers.index(kmer) == len(kmers) - 1:
                    records.append(record.id)
                    out_seq = str(record.seq[start:-(polya - 15)])
                    seqs.append(out_seq)
                    indices.append("NA")
                    found.append("NA")
                else:
                    continue
    elif type(trimmed_records) is str:
        with open(trimmed_records) as fasta:
            records_iter = SeqIO.parse(fasta, 'fasta')
            for record in records_iter:
                polya = int(lengths.get(record.id))
                start = len(record.seq) - (50 + polya)
                sliced_seq = str(record.seq[start:-polya])
                for kmer in kmers:
                    if kmer in sliced_seq:
                        index = sliced_seq.find(kmer)
                        index = index - 5
                        records.append(record.id)
                        out_seq = str(record.seq[start:-(polya - 15)])
                        seqs.append(out_seq)
                        indices.append(index)
                        found.append(kmer)
                        break
                    elif kmers.index(kmer) == len(kmers) - 1:
                        records.append(record.id)
                        out_seq = str(record.seq[start:-(polya - 15)])
                        seqs.append(out_seq)
                        indices.append("NA")
                        found.append("NA")
                    else:
                        continue


    df = pd.DataFrame(list(zip(records, seqs, found, indices)), columns=[
                      "Isoform", "3'_Sequence", "kmer", "Index"])
    output_path = pathlib.Path(output_directory)
    input_path = pathlib.Path(input_file)
    stem = input_path.stem
    with open(pathlib.Path.joinpath(output_path, stem + "_PAS.csv"), 'w') as pas:
        df.to_csv(pas)
    return df


def pssm(pas, output_directory, input_file):
    '''
    Retains a portion of the poly(A) tails, aligns by putative PAS, and
    generates a position-specific score matrix (PSSM) to calculate nucleotide
    frequencies at every position.
    '''
    sequences = {}

    if isinstance(pas, pd.DataFrame):
        data = pas
    elif type(pas) is str:
        data = pd.read_csv(pas)
    # create a dictionary from the read in CSV
    temp_seqs = {data["Isoform"].tolist()[i]: data["3'_Sequence"].tolist()[i]
                 for i in range(len(data["Isoform"].tolist()))}
    # create a dictionary where every sequence is the same length
    for id in temp_seqs:
        if type(temp_seqs[id]) is str and len(temp_seqs[id]) == 65:
            sequences[id] = temp_seqs[id]
    output_path = pathlib.Path(output_directory)
    if isinstance(pas, pd.DataFrame):
        input_path = pathlib.Path(input_file)
    elif type(pas) is str:
        input_path = pathlib.Path(pas)
    stem = input_path.stem
    # # write to FASTA the cleaned sequences
    with open(pathlib.Path.joinpath(output_path, stem + "_PAS_clean.fasta"), 'w') as output_file:
        for sequence in sequences:
            output_file.write(">" + sequence + "\n" +
                              sequences[sequence] + "\n")
    # read the FASTA back in
    alignment = AlignIO.read(pathlib.Path.joinpath(
        output_path, stem + "_pas_clean.fasta"), "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    pssm = summary_align.pos_specific_score_matrix(
        consensus, chars_to_ignore=["N"])
    sys.stdout = open(pathlib.Path.joinpath(
        output_path, stem + "_pssm.txt"), 'w')
    print(pssm)


def getargs():

    # arguments
    parser = argparse.ArgumentParser(
        description='Measure poly(A) tails, count k-mers, find the putative   \
        PAS, and generate PSSMs a given sequence file.')

    # required arguments
    parser.add_argument('input_file',
                        help='A FASTA file with long-read sequencing data    \
                        that hasn''t had poly(A) tails trimmed.')

    parser.add_argument('output_directory',
                        help='A path to the output directory.')

    # optional arguments
    parser.add_argument('-p', '--primer_file',
                        help='A file with a single primer sequences per      \
                        line. NOTE: it is recommended that the input         \
                        sequences are trimmed for adapters and primers prior \
                        to running this program. I do not guarantee a        \
                        desirable or predictable outcome with the -p         \
                        option.')

    parser.add_argument('-polyA', '--polyA', action="store_true",
                        help='Runs the algorithm to measure poly(A) tails in \
                        a FASTA file. Requires -p if sequences contain a 3\' \
                        adapter or primer. Will include poly(A) tail         \
                        measurement, k-mer counting in the 50 nt upstream of \
                        of the tail, identification of a putative PAS, and   \
                        generation of a PSSM for the entire input file.')

    parser.add_argument('-l', '--length_file',
                        help='A file with an isoform ID and poly(A) tail     \
                        length on each line to be used for k-mer counting    \
                        and PAS identification. This file can be generated   \
                        with -polyA.')

    parser.add_argument('-k', '--k', default=6, type=int,
                        help='An integer, the kmer length (defaults to 6).')

    parser.add_argument('-n', '--n', default=10, type=int,
                        help='An integer, the number of k-mers to pipe from  \
                        the k-mer search to the PAS search (defaults to 10).')

    parser.add_argument('-PAS', '--pas_file',
                        help='Identify the most likely PAS; requires a file  \
                        with one k-mer per line. NOTE: the algorithm will    \
                        search for each k-mer one at a time and stop the     \
                        search once it has found one. Thus, it is HIGHLY     \
                        recommended that this file is ordered in descending  \
                        PAS/k-mer usage.')

    parser.add_argument('-pssm', '--pssm',
                        help='Generate a position-specific score matrix for  \
                        the region flanking the PAS. Requires a CSV file     \
                        (generated by -polyA) with column names Isoform, 3\' \
                        Sequence, kmer, and Count.')

    args = parser.parse_args()
    return args


def main():

    args = getargs()

    if args.polyA:
        if args.primer_file:
            print("Primer file provided; trimming reads.")
            trim = primertrim(
                args.input_file, args.primer_file, args.output_directory)
            print("Measuring poly(A) tails.")
            lengths = measurepolya(
                trim, args.output_directory, args.input_file)
            print("Counting all k-mers in the 50 nt upstream of the start of the poly(A); retaining the %s most abundant." % args.n)
            kmer = findkmer(trim, lengths, args.output_directory,
                            args.k, args.n, args.input_file)
            print("Searching for the most likely PAS.")
            pas = findpas(trim, lengths, kmer,
                                 args.output_directory, args.input_file)
            print("Generating a position-specific score matrix for all transcripts.")
            pssm(pas, args.output_directory, args.input_file)
        elif not args.primer_file:
            print("Primer file not provided; assuming trimmed reads.")
            print("Measuring poly(A) tails.")
            lengths = measurepolya(
                args.input_file, args.output_directory,
                args.input_file)
            print("Counting all k-mers in the 50 nt upstream of the start of the poly(A); retaining the %s most abundant." % args.n)
            kmer = findkmer(args.input_file, lengths, args.output_directory,
                            args.k, args.n, args.input_file)
            print("Searching for the most likely PAS.")
            pas = findpas(args.input_file, lengths, kmer,
                                 args.output_directory, args.input_file)
            print("Generating a position-specific score matrix for all transcripts.")
            pssm(pas, args.output_directory, args.input_file)
    elif args.pssm:
        print("Generating a position-specific score matrix for provided transcripts.")
        pssm(args.pssm, args.output_directory, args.input_file)

main()
