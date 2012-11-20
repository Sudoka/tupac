#!/usr/bin/python -tt
import re
import sys
import cPickle as pickle

def main():
  c = pickle.load(open('data/chla_prot.pickle', 'r'))
  s = pickle.load(open('data/scen_gen.pickle', 'r'))

  count = 0
  for k in s:
    count += len(s[k])
  print 'Scenedesmus Genome'
  print '# of sequences =', len(s.keys())
  print '# of base pairs =', count, '\n'
  
  count = 0
  for k in c:
    for prot in c[k]:
      count += len(prot)
  print 'Chlamydomonas Proteome'
  print '# of sequences =', len(c.keys())
  print '# of amino acids =', count


# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()

