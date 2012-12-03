#!/usr/bin/python -tt
import re
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import cPickle as pickle

def main():
  
  seq_dict = pickle.load(open('data/chla_prot.pickle', 'r'))
  sample_seqs = [(k, seq_dict[k]) for k in seq_dict.keys()[:50]]
  s = write_OT(sample_seqs)
  open('sample.OT', 'w').write(s)

def write_OT(sample_seqs):
  s = ''
  query_string = ''
  for seq in sample_seqs:
    query_string += seq[0]+'\n'+seq[1]+'\n'

  blast_handle = NCBIWWW.qblast('tblastn', 'nr', query_string, entrez_query='scenedesmus dimorphus')
  blast_handle.seek(0)
  records = NCBIXML.parse(blast_handle)
  i = 0
  for record in records:
    if len(record.alignments) > 0:
      for align in record.alignments:
        row = sample_seqs[i][0] + '\t' + align.hit_id + '\t'
        frames = [hsp.frame[1] for hsp in align.hsps]
        if valid_align(frames):
          row += plus_or_minus(frames[0]) + '\t'
        else:
          row += '/' + '\t'
        query_coverage = float(sum([len(hsp.sbjct) for hsp in align.hsps])) / len(sample_seqs[i][1])
        if query_coverage < .9:
          row += str(1)
        else:
          row += str(2)
    else:
      row = sample_seqs[i][0] + '\t' + ' '*28 + '\t' + ' ' + '\t' + str(0)

    s += row + '\n'
    i += 1
  return s

def plus_or_minus(n):
  if n < 0:
    return '-'
  return '+'

#given list of frame values, returns true if all are < 0 or all are > 0
#this ensures the alignment is either sense or antisense derived
def valid_align(x):
  return all(item > 0 for item in x) or all(item < 0 for item in x)

# Standard boilerplate to call the main() function.
if __name__ == '__main__':
  main()

