class BoyerMoore:
    """Encapsulates pattern and its associated Boyer-Moore preprocessing"""
    
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.amap = {char: i for i, char in enumerate(self.alphabet)}
        # Make bad character rule table
        self.bad_char = self.dense_bad_char_tab(p, self.amap)
        # Make good suffix rule table
        self.good_suffix = self.good_suffix_table(p)
        self.full_shift = self.full_shift_table(p)
    
    def bad_character_rule(self, i, c):
        """Return # skips given by bad character rule at offset i"""
        assert c in self.amap
        ci = self.amap[c]
        assert i < len(self.bad_char)
        return self.bad_char[i][ci]
    
    def good_suffix_rule(self, i):
        """Return # skips given by good suffix rule at offset i"""
        return self.good_suffix[i]
    
    def match_skip(self):
        """Return # skips given by match skip rule"""
        return len(self.p) - self.full_shift[-1]
    
    def dense_bad_char_tab(self, p, amap):
        """Make dense bad character table"""
        tab = []
        nxt = [0] * len(amap)
        for i in range(len(p)):
            c = p[i]
            assert c in amap
            tab.append(list(nxt))
            nxt[amap[c]] = i + 1
        return tab
    
    def n_array(self, s):
        """Compute N array for string s"""
        n = [0] * len(s)
        for i in range(1, len(s)):
            k = i - 1
            while k >= 0 and s[k] == s[i + (i - 1 - k)]:
                k -= 1
            n[i] = i - k - 1
        return n
    
    def good_suffix_table(self, p):
        """Make good suffix table"""
        n = self.n_array(p)
        lp = [0] * len(p)
        for j in range(len(p) - 1):
            i = len(p) - n[j]
            if i < len(p):
                lp[i] = j + 1
        return lp
    
    def full_shift_table(self, p):
        """Make full shift table"""
        n = self.n_array(p)
        lp = [0] * len(p)
        for i in range(len(p)):
            lp[i] = n[i]
        for i in range(len(p) - 1, -1, -1):
            if lp[i] == 0:
                lp[i] = lp[i + 1] if i < len(p) - 1 else 0
        return lp

def boyer_moore_with_counts(p, p_bm, t):
    """
    Boyer-Moore algorithm that returns occurrences and statistics
    Returns: (occurrences, num_alignments, num_character_comparisons)
    """
    i = 0
    occurrences = []
    n_alignments = 0  # number of alignments tried
    n_comparisons = 0  # number of character comparisons
    
    while i < len(t) - len(p) + 1:
        n_alignments += 1
        shift = 1
        mismatched = False
        
        for j in range(len(p) - 1, -1, -1):
            n_comparisons += 1
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        
        i += shift
    
    return occurrences, n_alignments, n_comparisons

def readGenome(filename):
    """Parse FASTA file and return genome as a single string"""
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# Test with the provided example
def test_example():
    print("=" * 70)
    print("Testing with example:")
    print("=" * 70)
    p = 'word'
    t = 'there would have been a time for such a word'
    lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
    p_bm = BoyerMoore(p, lowercase_alphabet)
    occurrences, n_algn, n_char = boyer_moore_with_counts(p, p_bm, t)
    print(f"Pattern: '{p}'")
    print(f"Text: '{t}'")
    print(f"Result: {occurrences}, {n_algn}, {n_char}")
    print(f"Expected: [40], 12, 15")
    if occurrences == [40] and n_algn == 12 and n_char == 15:
        print("✓ TEST PASSED!")
    else:
        print("✗ TEST FAILED!")
    print()

# Main function to run on chr1.GRCh38.excerpt.fasta
def main():
    # Run test
    test_example()
    
    # Read the chromosome file
    print("=" * 70)
    print("Processing chr1.GRCh38.excerpt.fasta...")
    print("=" * 70)
    genome = readGenome('chr1.GRCh38.excerpt.fasta')
    print(f"Genome length: {len(genome)} bp")
    
    # The read to search for
    p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    print(f"Pattern length: {len(p)} bp")
    print(f"Pattern: {p}")
    print()
    
    # Create Boyer-Moore preprocessing
    p_bm = BoyerMoore(p, 'ACGT')
    
    # Run Boyer-Moore with counts
    occurrences, n_alignments, n_comparisons = boyer_moore_with_counts(p, p_bm, genome)
    
    print("Results:")
    print(f"  Occurrences found at positions: {occurrences}")
    print(f"  Number of alignments tried: {n_alignments}")
    print(f"  Number of character comparisons: {n_comparisons}")
    print("=" * 70)

if __name__ == "__main__":
    main()
