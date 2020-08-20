def FrequentWordsWithMismatchesAndReverseComplements(Text, k, d):
    FrequentPatterns = []
    Neighbourhoods = []
    Index = []
    Count = []
    l = len(Text)
    for i in range(l-k+1):
        Neighbourhoods.append(Neighbours(Text[i:i+k], d))
    NeighbourhoodArray = []
    for sublist in Neighbourhoods:
        for item in sublist:
            NeighbourhoodArray.append(item)
    NeighbourhoodArray2 = NeighbourhoodArray.copy()
    for string in NeighbourhoodArray2:
        NeighbourhoodArray.append(ReverseComplement(string))
    for i in range(len(NeighbourhoodArray)):
        Pattern = NeighbourhoodArray[i]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    Index.sort()
    for i in range(len(NeighbourhoodArray)-1):
        if Index[i] == Index[i+1]:
            Count[i+1] = Count[i] + 1
    maxCount = max(Count)
    for i in range(len(NeighbourhoodArray)):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(Index[i], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def FrequentWordsWithMismatches(Text, k, d):
    FrequentPatterns = []
    Neighbourhoods = []
    Index = []
    Count = []
    l = len(Text)
    for i in range(l-k+1):
        Neighbourhoods.append(Neighbours(Text[i:i+k], d))
    NeighbourhoodArray = []
    for sublist in Neighbourhoods:
        for item in sublist:
            NeighbourhoodArray.append(item)
    for i in range(len(NeighbourhoodArray)):
        Pattern = NeighbourhoodArray[i]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    Index.sort()
    for i in range(len(NeighbourhoodArray)-1):
        if Index[i] == Index[i+1]:
            Count[i+1] = Count[i] + 1
    maxCount = max(Count)
    for i in range(len(NeighbourhoodArray)):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(Index[i], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def ComputingFrequenciesWithMismatches(Text, k, d):
    FrequencyArray = []
    for i in range(4 ** k):
        FrequencyArray.append(0)
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i + k]
        Neighbourhood = Neighbours(Pattern, d)
        for ApproximatePattern in Neighbourhood:
            j = PatternToNumber(ApproximatePattern)
            FrequencyArray[j] += 1
    return FrequencyArray


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


def ImmediateNeighbours(Pattern):
    Neighbourhood = []
    Neighbourhood.append(Pattern)
    for i in range(len(Pattern)):
        symbol = Pattern[i]
        for x in 'ACGT':
            if x != symbol:
                Neighbour = Pattern[:i] + x + Pattern[i+1:]
                Neighbourhood.append(Neighbour)
    return Neighbourhood


def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i:i + len(Pattern)], Pattern) <= d:
            count = count + 1
    return count


def ApproximatePatternMatching(Pattern, Text, d):
    positions = []
    n = len(Text)
    k = len(Pattern)
    for i in range(n - k + 1):
        if HammingDistance(Pattern, Text[i:i + k]) <= d:
            positions.append(i)
    return positions


def HammingDistance(p, q):
    mismatch = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            mismatch += 1
    return mismatch


def MinimumSkew(Genome):
    positions = []
    skew = SkewArray(Genome)
    skew_min = min(skew.values())
    for key in skew:
        if skew[key] == skew_min:
            positions.append(key)
    return positions


def SkewArray(Genome):
    skew = {}
    Skew = []
    Skew.append(0)
    for i in range(len(Genome)):
        if Genome[i] == 'A' or Genome[i] == 'T':
            Skew.append(Skew[i])
        if Genome[i] == 'G':
            Skew.append(Skew[i] + 1)
        if Genome[i] == 'C':
            Skew.append(Skew[i] - 1)
    for i in range(len(Genome) + 1):
        skew[i] = Skew[i]
    return skew


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]

    array[0] = PatternCount(Genome[0:n // 2], symbol)

    for i in range(1, n):
        array[i] = array[i - 1]
        if ExtendedGenome[i - 1] == symbol:
            array[i] -= 1
        if ExtendedGenome[i + (n // 2) - 1] == symbol:
            array[i] += 1
    return array


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n // 2]
    for i in range(n):
        array[i] = PatternCount(ExtendedGenome[i:i + (n // 2)], symbol)
    return array


def PatternMatching(Pattern, Genome):
    positions = []
    n = len(Genome)
    k = len(Pattern)
    for i in range(n - k + 1):
        if Pattern == Genome[i:i + k]:
            positions.append(i)
    return positions


def ReverseComplement(Pattern):
    rev = ''
    rev = Reverse(Pattern)
    rev = Complement(rev)
    return rev


def Reverse(Pattern):
    rev = ''
    for char in Pattern:
        rev = char + rev
    return rev


def Complement(Pattern):
    rev = ''
    for char in Pattern:
        if char == 'A':
            rev = rev + 'T'
        elif char == 'T':
            rev = rev + 'A'
        elif char == 'G':
            rev = rev + 'C'
        elif char == 'C':
            rev = rev + 'G'
    return rev


def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    words.sort()
    return words


def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]
        freq[Pattern] = 0
    Pattern = Text[:k]
    for j in range(n - k + 1):
        Pattern = Text[j:j + k]
        freq[Pattern] += 1
    return freq


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count = count + 1
    return count


print(*FrequentWordsWithMismatchesAndReverseComplements('', 6, 3), sep=' ')
