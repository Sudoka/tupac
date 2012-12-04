#!/usr/bin/python -tt
import re
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import cPickle as pickle
import time

def main():
  
  seq_dict = pickle.load(open('data/chla_prot.pickle', 'r'))
  descs = sorted(pickle.load(open('data/three.descs', 'r')), key=lambda x: x[3], reverse=True)
  print len(descs)
  for desc in descs[:5]:
    print desc[3], '\n'
    print seq_dict[desc[0]], '\n'
    for align in desc[4]:
      for s in align[-1][:-1]:
        print s
      print ''
    print '\n\n'

  '''
  #offset for chunks of sequences to blast
  o = 0
  STEP_SIZE = 250
  start = time.time()
  seq_dict_keys = seq_dict.keys()
  n = len(seq_dict_keys)
  SECONDS_TO_RUN = 60*1
  while time.time()-start < SECONDS_TO_RUN and o+STEP_SIZE < n:
    sample_seqs = [(k, seq_dict[k]) for k in seq_dict_keys[o:o+STEP_SIZE]]
    descs = get_descriptions(sample_seqs)
    o += STEP_SIZE
    print time.time()-start
  pickle.dump(descs, open('data/one.descs', 'w'))
  '''

def get_descriptions(sample_seqs):
  query_string = ''
  for seq in sample_seqs:
    query_string += seq[0]+'\n'+seq[1]+'\n'

  blast_handle = NCBIWWW.qblast('tblastn', 'nr', query_string, entrez_query='scenedesmus dimorphus')
  blast_handle.seek(0)
  records = NCBIXML.parse(blast_handle)
  descs = []
  i = 0
  for record in records:
    if len(record.alignments) > 0:
      for align in record.alignments:
        desc = [sample_seqs[i][0], align.hit_id]

        frames = [hsp.frame[1] for hsp in align.hsps]
        if valid_align(frames):
          desc.append(plus_or_minus(frames[0]))
        else:
          desc.append('/')

        query_coverage = float(sum([len(hsp.sbjct) for hsp in align.hsps])) / len(sample_seqs[i][1])
        desc.append(query_coverage)

        #list of tuples of form ( (query_start, query_end), (sbjct_start, sbjct_end), (query, match, sbjct, frame) ) sorted by query_start
        hsp_info = sorted([((hsp.query_start, hsp.query_end), (hsp.sbjct_start, hsp.sbjct_end), (hsp.query, hsp.match, hsp.sbjct, hsp.frame[1])) for hsp in align.hsps], key= lambda t: t[0][0])
        desc.append(hsp_info)

    else:
      desc = [sample_seqs[i][0], ' ', ' ', 0.0, [], []]

    descs.append(desc)
    i += 1
  return descs

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

