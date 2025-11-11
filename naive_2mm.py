def readFastq(filename):
    """Parse FASTQ file and return sequences and quality strings"""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def readGenome(filename):
    """Parse FASTA file and return genome as a single string"""
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def phred33ToQ(qual):
    """Convert Phred+33 quality character to Q score"""
    return ord(qual) - 33

def naive_2mm(p, t):
    """
    Naïve matching algorithm allowing up to 2 mismatches
    Returns list of positions where pattern p matches text t with ≤2 mismatches
    """
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # mismatch
                mismatches += 1
                if mismatches > 2:  # more than 2 mismatches
                    break
        if mismatches <= 2:  # match with at most 2 mismatches
            occurrences.append(i)
    return occurrences

def main():
    # Read the reference genome
    genome = readGenome('lambda_virus.fa')
    
    # Read the FASTQ file
    sequences, qualities = readFastq('myReads.fastq')
    
    # Counters
    reads_passed = 0
    reads_failed = 0
    
    print("Processing reads...")
    print("-" * 60)
    
    for i, (seq, qual) in enumerate(zip(sequences, qualities)):
        # Check if all base qualities are > 5
        quality_scores = [phred33ToQ(q) for q in qual]
        
        if all(q > 5 for q in quality_scores):
            reads_passed += 1
            # Find matches with up to 2 mismatches
            matches = naive_2mm(seq, genome)
            
            if matches:
                print(f"Read {i}: Found at position(s) {matches}")
            else:
                print(f"Read {i}: No match found")
        else:
            reads_failed += 1
    
    print("-" * 60)
    print(f"\nSummary:")
    print(f"Reads that passed quality filter (all bases Q > 5): {reads_passed}")
    print(f"Reads that failed quality filter: {reads_failed}")
    print(f"Total reads: {reads_passed + reads_failed}")

if __name__ == "__main__":
    main()
