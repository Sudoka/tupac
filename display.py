#!/usr/bin/python -tt
import re
from sys import argv as args
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import cPickle as pickle
import time

def main():
  # Arguments:
  # -sf    -  pickled seq_dict file   
  # -df    -  pickled descriptor file
  # -info  -  what type of info you want, options include:
  #             -orthotable
  #             -exontable
  #             -general
  # -of    -  outfile to write the info to, if supplied with 'stdout' it will print to screen
  # -t     -  integer argument specifying top hits to write, if -1 will write all hits in unsorted order

  if '-h' in args or '--h' in args or '--help' in args:
    print 'Usage is as follows: ./display.py -sf seq_dict_file -df descriptor_file -info info_type -of outfile -t top_hits'
    sys.exit()
  vals = [args[args.index(flag) + 1] for flag in ['-sf', '-df', '-info', '-of', '-t']]

  info = vals[2]
  outfile = vals[3]
  top_hits = int(vals[4])
  seq_dict = pickle.load(open(vals[0], 'r'))
  descs = pickle.load(open(vals[1], 'r'))

  if top_hits != -1:
    descs = sorted(descs, key= lambda x: x[3], reverse=True)
    descs = descs[:top_hits]

  print_str = ''

  if info == 'orthotable':
    for desc in descs:
      print_str += orthotable_row(desc)
  elif info == 'exontable':
    for desc in descs:
      print_str += exontable_row(desc)
  elif info == 'general':
    for desc in descs:
      print_str += general_row(desc)

  if outfile == 'stdout':
    print print_str
  else:
    open(outfile, 'w').write(print_str)


def orthotable_row(desc):
  if desc[3] == 0.0:
    return desc[0] + '\t' + 'NO_HIT_FOUND' + ' '*16 +'\t \t' + str(0) + '\n'
  if desc[3] > .9:
    qc = 2
  else:
    qc = 1
  return desc[0] + '\t' + desc[1] + '\t' + desc[2] + '\t' + str(qc) + '\n'

def exontable_row(desc):
  if desc[3] == 0.0:
    return desc[0] + '\t' + 'NO_HITS_FOUND'+' '*16 + '\n'
  return desc[0] + '\t' + desc[1] + '\t' + ''.join( [ str(align[1]) + ' ' for align in desc[4] ] ) + '\n'

def general_row(desc):
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

