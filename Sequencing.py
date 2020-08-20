import random

def StringSpelledByGappedPatterns(GappedPatterns, k, d):
    FirstPatterns = []
    SecondPatterns = []
    for pattern in GappedPatterns:
        FirstPatterns.append(pattern[0:k])
    for pattern in GappedPatterns:
        SecondPatterns.append(pattern[k+1:2*k+1])
    PrefixString = PathToGenome(FirstPatterns)
    SuffixString = PathToGenome(SecondPatterns)
    print(FirstPatterns)
#    print(PrefixString)
    print(SecondPatterns)
#    print(SuffixString)
    l = len(PrefixString)
    for i in range(k+d+1, l):
        if PrefixString[i] != SuffixString[i-k-d]:
            return 'there is no string spelled by the gapped patterns'
    return PrefixString + SuffixString[-k-d:]


def PairedComposition(k, d, Text):
    pairs = []
    l = len(Text)
    for i in range(l-k-d-k+1):
        kmer = Text[i:i+k] + '|' + Text[i+k+d:i+k+d+k]
        pairs.append(kmer)
    pairs.sort()
    return pairs


def kUniversalCircularString(k):
    strings = BinaryStrings(k)
    String_with_edges = StringReconstruction(strings)
    circular = String_with_edges[:-(k-1)]
    return circular


def BinaryStrings(k):
    strings = []
    if k == 0:
        return strings
    elif k > 0:
        strings = []
        for i in range(2**k):
            string = bin(i).replace("0b","")
            while len(string) != k:
                string = '0' + string
            strings.append(string)
        return strings


def StringReconstruction(Patterns):
    dB = DeBruijin2(Patterns)
    path = EulerianPath(dB)
    #path = EulerianCycle(dB)
    Text = PathToGenome(path)
    return Text


def EulerianPath(Graph):
    stack = []
    circuit = []
    list_keys = []
    list_values = []
    for key in Graph:
        list_keys.append(key)
        list_values.append(Graph[key])
    flat_list = []
    for sublist in list_values:
        for item in sublist:
            flat_list.append(item)
    list_values = flat_list
    missing = ''
    for item in list_values:
        if item not in list_keys:
            missing = item
    outdeg = {}
    for key in Graph:
        outdeg[key] = len(Graph[key])
    outdeg[missing] = 0
    missing2 = ''
    for item in list_keys:
        if item not in list_values:
            missing2 = item
    indeg = {}
    for item in list_values:
        indeg[item] = 0
    for item in list_values:
        indeg[item] += 1
    indeg[missing2] = 0
    full_list = list_values + list_keys
    full_list = list(dict.fromkeys(full_list))
    for item in full_list:
        if outdeg[item] - indeg[item] == 1:
            beginning = item
        elif indeg[item] - outdeg[item] == 1:
            ending = item
    if ending not in Graph:
        Graph.update({ending:[beginning]})
    elif ending in Graph:
        Graph[ending].append(beginning)
    location = beginning
    while len(stack) != 0 or len(Graph[location]) != 0:
        if len(Graph[location]) == 0:
            circuit.append(location)
            location = stack[-1]
            del stack[-1]
        else:
            stack.append(location)
            new_location = random.choice(Graph[location])
            Graph[location].remove(new_location)
            location = new_location
    circuit.reverse()
    while circuit[0] != beginning or circuit[-1] != ending:
        circuit.insert(0, circuit[-1])
        del circuit[-1]
    return circuit


def EulerianCycle(Graph):
    stack = []
    circuit = []
    location = random.choice(list(Graph))
    while len(stack) != 0 or len(Graph[location]) != 0:
        if len(Graph[location]) == 0:
            circuit.append(location)
            location = stack[-1]
            del stack[-1]
        else:
            stack.append(location)
            new_location = random.choice(Graph[location])
            Graph[location].remove(new_location)
            location = new_location
    circuit.reverse()
    circuit.append(circuit[0])
    return circuit


def DeBruijin2(Patterns):
    debruijin = {}
    l = len(Patterns)
    k = len(Patterns[0])
    for kmer in Patterns:
        prefix = Prefix(kmer)
        suffix = Suffix(kmer)
        if prefix not in debruijin:
            debruijin.update({prefix:[suffix]})
        elif prefix in debruijin:
            debruijin[prefix].append(suffix)
    return debruijin


def DeBruijin(k, Text):
    debruijin = {}
    kmers = Composition(k-1, Text)
    d = len(kmers)
    for i in range(d-1):
        if kmers[i] not in debruijin:
            debruijin.update({kmers[i]:[kmers[i+1]]})
        elif kmers[i] in debruijin:
            debruijin[kmers[i]].append(kmers[i+1])
    return debruijin


def Overlap(Patterns):
    overlap = {}
    for i in Patterns:
        for j in Patterns:
            if Suffix(i) == Prefix(j):
                if i not in overlap:
                    overlap.update({i:[j]})
                elif i in overlap:
                    overlap[i].append(j)
    return overlap


def Prefix(Pattern):
    return Pattern[0:len(Pattern)-1]


def Suffix(Pattern):
    return Pattern[1:len(Pattern)]


def PathToGenome(Path):
    Genome = ''
    size = len(Path[0])
    i = len(Path)
    Genome += Path[0]
    for j in range(i-1):
        Genome += Path[j+1][size-1]
    return Genome

def Composition(k, Text):
    kmers = []
    d = len(Text)
    for i in range(d-k+1):
        kmer = Text[i:i+k]
        kmers.append(kmer)
    #comment out for DeBruijin(k, Text):
    #random.shuffle(kmers)
    return kmers

with open('names.txt', 'r') as f:
    GappedPatterns = [line.strip() for line in f]

#code for Overlap(Patterns)
#overlap = Overlap(Patterns)
#for i in overlap:
#    print(i + ' -> ', end = '')
#    print(*overlap[i], sep = ',')


#code for DeBruijin(k, Text)/DeBruijin2(Patterns)
#debruijin = DeBruijin2(Patterns)
#for i in debruijin:
#    print(i + ' -> ', end = '')
#    print(*debruijin[i], sep = ',')


#code for EulerianCycle(Graph)/EulerianPath(Graph)
#with open('names.txt', 'r') as file:
#    Graph = dict((line.strip().split(' -> ') for line in file))
#    for key in Graph:
#        Graph[key] = Graph[key].split(',')

print(DeBruijin2(GappedPatterns))