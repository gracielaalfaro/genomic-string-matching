class BoyerMoore:
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        self.amap = {c: i for i, c in enumerate(alphabet)}
        self.bad_char = self.bad_character_table(p, self.amap)
        self.n = self.n_array(p)
        self.l_prime = self.good_suffix_table(p, self.n)
        self.l_prime_rev = self.match_skip_table(p, self.n)

    # ---------- BAD CHARACTER RULE ----------
    def bad_character_table(self, p, amap):
        n = len(amap)
        m = len(p)
        occ = [[-1] * m for _ in range(n)]
        for i in range(m):
            c = p[i]
            if c in amap:
                ci = amap[c]
                for j in range(i, m):
                    occ[ci][j] = i
        return occ

    def bad_character_rule(self, i, c):
        if c not in self.amap:
            return i + 1
        ci = self.amap[c]
        if self.bad_char[ci][i] == -1:
            return i + 1
        return i - self.bad_char[ci][i]

    # ---------- GOOD SUFFIX RULE ----------
    def z_array(self, s):
        Z = [0] * len(s)
        Z[0] = len(s)
        right, left = 0, 0
        for k in range(1, len(s)):
            if k > right:
                n = 0
                while n + k < len(s) and s[n] == s[n + k]:
                    n += 1
                Z[k] = n
                if n > 0:
                    left = k
                    right = k + n - 1
            else:
                k1 = k - left
                beta_len = right - k + 1
                if Z[k1] < beta_len:
                    Z[k] = Z[k1]
                else:
                    n = right + 1
                    while n < len(s) and s[n] == s[n - k]:
                        n += 1
                    Z[k] = n - k
                    left = k
                    right = n - 1
        return Z

    def n_array(self, s):
        """Generate N array (length of longest suffix of s[i:] matching a prefix of s)."""
        rev_s = s[::-1]
        z_rev = self.z_array(rev_s)
        return z_rev[::-1]

    def good_suffix_table(self, p, n):
        m = len(p)
        l_prime = [0] * m
        for j in range(m - 1):
            i = m - n[j]
            if i < m:
                l_prime[i] = j + 1
        return l_prime

    def match_skip_table(self, p, n):
        m = len(p)
        lp = [0] * (m + 1)
        longest = 0
        for i in range(m - 1, -1, -1):
            if n[i] == i + 1:
                longest = n[i]
            lp[m - i - 1] = longest
        return lp

    def good_suffix_rule(self, j):
        m = len(self.p)
        if j == m - 1:
            return 0
        i = j + 1
        if self.l_prime[i] > 0:
            return m - self.l_prime[i] - 1
        return m - self.l_prime_rev[i]

    def match_skip(self):
        return len(self.p) - self.l_prime_rev[1]

# ---------- MAIN FUNCTION ----------
def boyer_moore_with_counts(p, p_bm, t):
    occurrences = []
    n_alignments = 0
    n_comparisons = 0
    i = 0
    while i <= len(t) - len(p):
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
            shift = max(shift, p_bm.match_skip())
        i += shift
    return occurrences, n_alignments, n_comparisons

# ---------- TEST HARNESS ----------
def main():
    print("=" * 70)
    print("Testing with example:")
    print("=" * 70)
    p = "word"
    t = "there would have been a time for such a word"
    alphabet = "abcdefghijklmnopqrstuvwxyz "
    p_bm = BoyerMoore(p, alphabet)
    occurrences, n_align, n_char = boyer_moore_with_counts(p, p_bm, t)
    print(f"Pattern: '{p}'")
    print(f"Text: '{t}'")
    print(f"Result: {occurrences}, {n_align}, {n_char}")
    print("Expected: [40], 12, 15")
    if occurrences == [40] and n_align == 12 and n_char == 15:
        print("✓ TEST PASSED!")
    else:
        print("✗ TEST FAILED!")
    print()

    # Genome test (optional)
    print("=" * 70)
    print("Processing chr1.GRCh38.excerpt.fasta...")
    print("=" * 70)
    try:
        with open("chr1.GRCh38.excerpt.fasta") as f:
            lines = f.readlines()
        genome = ''.join([line.strip() for line in lines if not line.startswith('>')])
        print(f"Genome length: {len(genome)} bp")

        p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
        print(f"Pattern length: {len(p)} bp")
        print(f"Pattern: {p}")

        p_bm = BoyerMoore(p, "ACGT")
        occurrences, n_align, n_char = boyer_moore_with_counts(p, p_bm, genome)
        print(f"\nOccurrences: {occurrences}")
        print(f"Alignments tried: {n_align}")
        print(f"Character comparisons: {n_char}")
    except FileNotFoundError:
        print("chr1.GRCh38.excerpt.fasta not found — skipping genome test.")

if __name__ == "__main__":
    main()
