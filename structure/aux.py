#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys

def levicivita(a,b,c):
  '''
  levi civity tensor hard coded for python input in range [0,1,2]
  '''
  a+=1
  b+=1
  c+=1
  l = [a,b,c]
  s = set(l)
  if min(l) < 1 or max(l) > 3: return 0
  if len(s) != 3: return 0
  if a==1:
    if b==2: return 1  # 123
    else: return -1    # 132
  if a==2:
    if b==1: return -1 # 213
    else: return 1     # 231
  if a==3:
    if b==1: return 1  # 312
    else: return -1    # 321

# The MIT License (MIT)
# Copyright (c) 2016 Vladimir Ignatev
# https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
def progressBar(count, total, status='', prefix=''):
  bar_len = 60
  filled_len = int(round(bar_len * count / float(total)))
  percents = round(100.0 * count / float(total), 1)
  bar = '=' * filled_len + '-' * (bar_len - filled_len)
  sys.stdout.write('{} [{}] {}{} ... {}\r'.format(prefix, bar, percents, '%', status))
  sys.stdout.flush()
  if count == total:
    sys.stdout.write('\n')
    sys.stdout.flush()


if __name__ == '__main__':
  for i in range(3):
    for j in range(3):
      for k in range(3):
        print('e_',i,j,k,' = ', levicivita(i,j,k))
