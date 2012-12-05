#!/usr/bin/python -tt
import re
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import cPickle as pickle
import time

def main():
  
  seq_dict = pickle.load(open('data/chla_prot.pickle', 'r'))
  descs = sorted(pickle.load(open('data/20minutes.desc', 'r')), key=lambda x: x[3], reverse=True)
  for desc in descs[:5]:
    print homolog_info(desc)

def homolog_info(desc):
  s = 'Protein ID:\t'+desc[0]+'\n'
  s += 'Gene ID:\t'+desc[1]+'\n'
  s += 'Query Coverage:\t'+str(int(desc[3]*100))+'%'+'\n'
  s += 'Senses:\t'+str([p_or_m(align[2][3]) for align in desc[4]])+'\n'
  s += 'Query / Gene Coords:\t'
  for align in desc[4]:
    #s += str(align[0][0])+':'+str(align[0][1])+'\t'
    s += 'Q: ' + str(align[0]) + ' G: ' + str(align[1]) + '\t'
  s += '\n\n'
  s += pretty_print(align_strings(desc))
  return s

def pretty_print(align_strings):
  n = len(align_strings[0])
  offset = 0
  s = ''
  while offset+80 < n:
    s += align_strings[0][offset:offset+80]+'\n'
    s += align_strings[1][offset:offset+80]+'\n'
    s += align_strings[2][offset:offset+80]+'\n\n'
    offset += 80
  if not align_strings[0][offset:] == '':
    s += align_strings[0][offset:]+'\n'
    s += align_strings[1][offset:]+'\n'
    s += align_strings[2][offset:]+'\n\n'
  return s

def align_strings(desc):
  s1 = ''; s2 = ''; s3 = '';
  for aligns in desc[4]:
    s1 += '         ' + aligns[2][0]
    s2 += '   ~~~   ' + aligns[2][1]
    s3 += '         ' + aligns[2][2]
  return (s1[9:], s2[9:], s3[9:])
      
def p_or_m(n):
  if n < 0:
    return '-'
  return '+'


# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()

