import sys
from collections import defaultdict
from operator import itemgetter

bedfile = sys.argv[1]
output = sys.argv[2]

outdictionary = defaultdict(list)
with open(bedfile,'r') as bedfileopen:
  for line in bedfileopen:
    linesplit = line.strip('\n').split('\t')
    ident = linesplit[3]
    strand = linesplit[5]
    start = linesplit[1]
    if len(linesplit) == 18:
      gene = linesplit[15]
    elif len(linesplit) == 21:
      gene = linesplit[20]
    outdictionary[ident].append((int(start),strand,gene,line))


dropped = 0
spliced = 0
with open(output,'w') as outputopen:     
  for k,v in outdictionary.iteritems():
    if len(v) == 1:
      outputopen.write(str(v[0][3]))
    elif len(v) == 2:
      if v[0][2] == v[1][2]:
        spliced += 1
        v.sort(key = itemgetter(0))
        if v[0][1] == '+':
          outputopen.write(str(v[0][3]))
        elif v[0][1] == '-':
          outputopen.write(str(v[1][3]))
      else:
        dropped += len(v)
    elif len(v) > 2:
      dropped += len(v)

print 'dropped: '+str(dropped)
print 'spliced: '+str(spliced)
