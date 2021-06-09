#to filter intersect with proportion less then 0.5
from __future__ import division
import sys
for line in open(sys.argv[1],"r") :
        line = line.strip()
        info = line.split("\t")
        a_len = sum([int(i) for i in info[10].split(",")])
        b_len = sum([int(i) for i in info[22].split(",")])
        o_len = int(info[24])
        if (int(info[1])-int(info[13]))*(int(info[2])-int(info[14])) <= 0 :
                print line
        else:
                if (o_len/a_len) >= 0.5 or (o_len/b_len) >= 0.5 :
                        print line
