import random

def MedianString(Dna, k):
    distance = 2147483646
    for i in range(4**k):
        Pattern = NumberToPattern(i, k)
        if distance > DistanceBetweenPattenAndStrings(Pattern, Dna):
            distance = DistanceBetweenPattenAndStrings(Pattern, Dna)
            Median = Pattern
    return Median


def NumberToPattern(index, k):
    if k == 1:
        return NumberToSymbol(index)
    prefixIndex = Quotient(index, 4)
    r = Reminder(index, 4)
    symbol = NumberToSymbol(r)
    PrefixPattern = NumberToPattern(prefixIndex, k - 1)
    return PrefixPattern + symbol


def NumberToSymbol(Number):
    if Number == 0:
        return "A"
    elif Number == 1:
        return "C"
    elif Number == 2:
        return "G"
    elif Number == 3:
        return "T"


def Reminder(n, m):
    reminder = n % m
    return reminder


def Quotient(n, m):
    quotient = n // m
    return quotient


def PatternToNumber(Pattern):
    if Pattern == "":
        return 0
    symbol = LastSymbol(Pattern)
    prefix = Prefix(Pattern)
    number = 4 * PatternToNumber(prefix) + SymbolToNumber(symbol)
    return number


def SymbolToNumber(Symbol):
    if Symbol == "A":
        return 0
    elif Symbol == "C":
        return 1
    elif Symbol == "G":
        return 2
    elif Symbol == "T":
        return 3


def Prefix(Pattern):
    prefix = ""
    t = len(Pattern)
    prefix = Pattern[0:t - 1]
    return prefix


def LastSymbol(Pattern):
    symbol = ""
    t = len(Pattern)
    symbol = Pattern[t - 1]
    return symbol


def DistanceBetweenPattenAndStrings(Pattern, Dna):
    k = len(Pattern)
    distance = 0
    for string in Dna:
        Hammingdistance = 2147483646
        kmers = []
        for i in range(len(string) - k + 1):
            kmers.append(string[i:i + k])
        for kmer in kmers:
            if Hammingdistance > HammingDistance(Pattern, kmer):
                Hammingdistance = HammingDistance(Pattern, kmer)
        distance += Hammingdistance
    return distance


def MotifEnumeration(Dna, k, d):
    Patterns = []
    World = []
    for string in Dna:
        kmers = []
        for i in range(len(string)-k+1):
            kmers.append(string[i:i+k])
        neighbourhood = set()
        for kmer in kmers:
            neighbours = Neighbours(kmer, d)
            for kmer2 in neighbours:
                neighbourhood.add(kmer2)
        neighbourhood_list = list(neighbourhood)
        World.append(neighbourhood_list)
    for i in range(len(World)):
        for j in range(len(World[i])):
            # change depending on the len(Dna)
            if (World[i][j] in World[0]) and (World[i][j] in World[1]) and (World[i][j] in World[2]) and (World[i][j] in World[3]):
                Patterns.append(World[i][j])
    Patterns2 = set(Patterns)
    return Patterns2


def Neighbours(Pattern, d):
    if d == 0:
        return [Pattern]
    if len(Pattern) == 1:
        return ['A', 'C', 'G', 'T']
    Neighbourhood = []
    SuffixNeighbours = Neighbours(Suffix(Pattern), d)
    for Text in SuffixNeighbours:
        if HammingDistance(Suffix(Pattern), Text) < d:
            for x in 'ACGT':
                Neighbourhood.append(x + Text)
        else:
            Neighbourhood.append(FirstSymbol(Pattern) + Text)
    return Neighbourhood


def FirstSymbol(Pattern):
    first = Pattern[0]
    return first


def Suffix(Pattern):
    suffix = Pattern[1:]
    return suffix


def HammingDistance(p, q):
    mismatch = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            mismatch += 1
    return mismatch


# Input:  Integers k, t, and N, followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
# !!!!! Change in def Consensus(Motifs):
# the line:
#       count = Count(Motifs)
# to the line
#       count = CountWithPseudocounts(Motifs)
def GibbsSampler(Dna, k, t, N):
    BestMotifs = [] # output variable
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(N):
        i = random.randint(0, t-1)
        Motifs.pop(i)
        Profile = ProfileWithPseudocounts(Motifs)
        Motifs.insert(i, ProfileGeneratedString(Dna[i], Profile, k))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {}
    for i in range(0, n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = '' # output variable
    n = random.uniform(0,1)
    for kmer in Probabilities:
        n -= Probabilities[kmer]
        if n <= 0:
            return kmer


# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    normalized = {}
    sum = 0
    for key in Probabilities:
        sum += Probabilities[key]
    for key in Probabilities:
        normalized[key] = Probabilities[key] / sum
    return normalized


# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs
        

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    kmers = []
    s = len(Dna[0])
    for string in Dna:
        m = random.randint(0, s - k)
        kmer = string[m:m + k]
        kmers.append(kmer)
    return kmers


# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    kmers = []
    k = len(Profile['A'])
    for string in Dna:
        kmer = ProfileMostProbableKmer(string, k, Profile)
        kmers.append(kmer)
    return kmers


# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
# !!!!! Change in def Consensus(Motifs):
# the line:
#       count = Count(Motifs)
# to the line
#       count = CountWithPseudocounts(Motifs)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    k = len(Motifs[0])
    profile = {}
    profile = CountWithPseudocounts(Motifs)
    for key in profile:
        for j in range(k):
            profile[key][j] = profile[key][j] / t
    return profile


# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {}  # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


def ProfileMostProbableKmer(text, k, profile):
    kmer = text[0:k]
    t = len(text)
    m = 0
    for i in range(t - k + 1):
        pr = Pr(text[i:i + k], profile)
        if pr > m:
            m = pr
            kmer = text[i:i + k]
    return kmer


# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    p = 1
    t = len(Text)
    for i in range(t):
        for key in Profile:
            if Text[i] == key:
                p = p * Profile[Text[i]][i]
    return p


# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    score = 0
    consensus = Consensus(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score


# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    #count = Count(Motifs)
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    for key in profile:
        for j in range(k):
            profile[key][j] = profile[key][j] / t
    return profile


# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {}  # initializing the count dictionary
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

with open('names.txt', 'r') as f:
    Dna = [line.strip() for line in f]

t = len(Dna)
k = 15
N = 1500
# Call GibbsSampler(Dna, k, t, N) 20 times and store the best output in a variable called BestMotifs
BestMotifs = GibbsSampler(Dna, k, t, N)
for i in range(19):
    result = GibbsSampler(Dna, k, t, N)
    if Score(result) < Score(BestMotifs):
        BestMotifs = result

# Print the BestMotifs variable
print(*BestMotifs, sep='\n')
