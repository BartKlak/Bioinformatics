def BetterClumpFinding(Genome, k, L, t):
    FrequentPatterns = []
    Clump = []
    d = len(Genome)
    for i in range(4**k):
        Clump.append(0)
    Text = Genome[0:L]
    FrequencyArray = ComputingFrequencies(Text, k)
    for index in range(4 ** k):
        if FrequencyArray[index] >= t:
            Clump[index] = 1
    for i in range(1, d-L+1):
        FirstPattern = Genome[i-1:i-1+k]
        index = PatternToNumber(FirstPattern)
        FrequencyArray[index] -= 1
        LastPattern = Genome[i+L-k: i+L]
        index = PatternToNumber(LastPattern)
        FrequencyArray[index] += 1
        if FrequencyArray[index] >= t:
            Clump[index] = 1
    for i in range(4**k):
        if Clump[i] == 1:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def ClumpFinding(Genome, k, L, t):
    FrequentPatterns = []
    Clump = []
    d = len(Genome)
    for i in range(4**k):
        Clump.append(0)
    for i in range(d-L+1):
        Text = Genome[i:i+L]
        FrequencyArray = ComputingFrequencies(Text, k)
        for index in range(4**k):
            if FrequencyArray[index] >= t:
                Clump[index] = 1
    for i in range(4**k):
        if Clump[i] == 1:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def FindingFrequentWordsBySorting(Text, k):
    FrequentPatterns = []
    Index = []
    Count = []
    l = len(Text)
    for i in range(l-k+1):
        Pattern = Text[i:i+k]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    Index.sort()
    for i in range(1, l-k+1):
        if Index[i] == Index[i-1]:
            Count[i] = Count[i-1] + 1
    maxCount = max(Count)
    for i in range(l-k+1):
        if Count[i] == maxCount:
            Pattern = NumberToPattern(Index[i], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def FasterFrequentWords(Text, k):
    FrequentPatterns = []
    FrequencyArray = ComputingFrequencies(Text, k)
    maxCount = max(FrequencyArray)
    for i in range(4 ** k):
        if FrequencyArray[i] == maxCount:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns


def ComputingFrequencies(Text, k):
    FrequencyArray = []
    for i in range(4 ** k):
        FrequencyArray.append(0)
    for i in range(len(Text) - k + 1):
        Pattern = Text[i:i + k]
        j = PatternToNumber(Pattern)
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



print(FindingFrequentWordsBySorting('ATATATATATAGAGCACACCCCAACACCA', 3))