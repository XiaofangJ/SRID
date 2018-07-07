#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
import tempfile

def run_cmd(cmd):
    cmd = "t1=`date +%s`;"+cmd+";t2=`date +%s`;tdiff=`echo 'scale=3;('$t2'-'$t1')/60' | bc`;echo '##### Total time:  '$tdiff' mins'"
    p = subprocess.Popen(cmd,bufsize=-1, shell=True, universal_newlines=True, stdout=subprocess.PIPE,executable='/bin/bash')
    output = p.communicate()[0]
    return output

def is_tool(name):
    """Check whether `name` is on PATH."""
    from distutils.spawn import find_executable
    return find_executable(name) is not None

def pipeline(args):
    input=args.input
    output=args.output
    core=args.core
    min_SR=args.min_SR
    mean=args.mean
    sd=args.sd
    readlen=args.readlen
    tmpdir=args.tmpdir

    f=next(tempfile._get_candidate_names())
    prefix=tmpdir+f

    # step 1: identify read pairs that meet the following requirements:
    #        a. Both reads aligned but not in proper pair (-F 14)
    #        b. First in pair (Flag == 81 or 97)
    #        c. Both reads align to the same scaffold ( $7== "=")
    #        d. The entire read of the first in pair can aligned ($6 ~/^[0-9]*M$/)
    cmd='''
    samtools view -F 14 -@ {core} {input}|awk 'BEGIN{{OFS="\\t"}}{{if(($2==81||$2==97) && $7 == "=" && $6~/^[0-9]*M$/){{if($4<$8){{print $3,$4,$8,$1}}else{{print $3,$8,$4,$1}}}}}}'|sortBed >{prefix}.PE.bed '''.format(input=input,prefix=prefix,core=core)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)

    # step 2: identify reads that aligned to two position of the reference genomes
    #      2.1 Remove reads that don't aligned to the reference genomes and sort the sam file based on reads name
    #      2.2 Retrieve reads that have two entries in the sam files and store those entries
    cmd = '''
    samtools view -F 4 -hb -@ {core} {input} |samtools sort -n -@ {core} > {prefix}.tmp.bam
    samtools view -f 65 -@ {core} {prefix}.tmp.bam|sort -k1,1 >{prefix}.tmp.1.sam;
    cut -f 1 {prefix}.tmp.1.sam |uniq -dc|awk '$1==2{{print $2}}'|sort -k1,1 >{prefix}.1.duplist;
    join -t $'\\t' {prefix}.1.duplist  {prefix}.tmp.1.sam >{prefix}.1.sam 
    
    samtools view -f 129  -@ {core} {prefix}.tmp.bam|sort -k1,1 >{prefix}.tmp.2.sam;
    cut -f 1 {prefix}.tmp.2.sam |uniq -dc|awk '$1==2{{print $2}}'|sort -k1,1 >{prefix}.2.duplist;
    join -t $'\\t' {prefix}.2.duplist  {prefix}.tmp.2.sam >{prefix}.2.sam '''.format(input=input,core=core, min_SR=min_SR,prefix=prefix)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)

    # step 3: identify split reads and
    #        3.1. Convert sam to bed format, and sort bed file based on reads name and coordinates where the reads mapped in the genomes
    #        3.2. Reads are split reads that supporting insertions in reference genomes:
    #           a. Both split fragments aligned to the same scaffold
    #           b. The first split fragment matched in pattern ^([0-9]*)M([0-9]*)[SH]$, while the second fragment match in pattern match($19,/^([0-9]*)[SH]([0-9]*)M$/
    #           c. There is no base pair that failed to align to the scaffold ($4>=0)
    #           d. More than min_SR number of reads supporting the split site
    #           e. The insertions region is with range 1000 ~ 150000 bps
    cmd = '''
    sam2bed <{prefix}.1.sam |sort -k4,4 -k2,2n >{prefix}.1.bed;
    awk 'BEGIN{{OFS="\\t";ORS=""}}NR%2==1{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"\\t"}}NR%2==0{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"\\n"}}' {prefix}.1.bed |
    awk '$1==$12 && match($8,/^([0-9]*)M([0-9]*)[SH]$/,m) && match($19,/^([0-9]*)[SH]([0-9]*)M$/,n) && m[1] >= 10 && n[2] >= 10 && m[1]-n[1]>0 {{print $1 "\\t" $2+m[0] "\\t" $13 "\\t" $4}}'|sed 's/$/.1/' > {prefix}.1.SR.bed

    sam2bed <{prefix}.2.sam |sort -k4,4 -k2,2n >{prefix}.2.bed;
    awk 'BEGIN{{OFS="\\t";ORS=""}}NR%2==1{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"\\t"}}NR%2==0{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11"\\n"}}' {prefix}.2.bed |
    awk '$1==$12 && match($8,/^([0-9]*)M([0-9]*)[SH]$/,m) && match($19,/^([0-9]*)[SH]([0-9]*)M$/,n) && m[1] >= 10 && n[2] >= 10 && m[1]-n[1]>0 {{print $1 "\\t" $2+m[0] "\\t" $13 "\\t" $4}}'|sed 's/$/.2/' > {prefix}.2.SR.bed

    cat {prefix}.1.SR.bed {prefix}.2.SR.bed >{prefix}.SR.bed

    cut -f 1-3 {prefix}.SR.bed |sort|uniq -c|awk '$1>={min_SR}{{print $2"\\t"$3"\\t"$4}}'| awk '$3-$2>1000 && $3-$2<150000' |sortBed >{prefix}.splitsite.bed '''.format(input=input,core=core, min_SR=min_SR,prefix=prefix)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)

    # step 4: identify reads pairs that support the split site identify in step 3 based on insertion size
    # the insertion size - insertions size should be with the range of the insertion size estimated from reads aligned proper
    cmd= '''
    closestBed -a {prefix}.splitsite.bed -b {prefix}.PE.bed|
    awk 'BEGIN{{OFS="\\t"}}{{insize=$6-$5+{readlen}-$3+$2;if(insize >= {mean} - 2*{sd} && insize <= {mean}+2*{sd}){{print $1,$2,$3,$7}}}}' > {prefix}.support.PE.bed '''.format(prefix=prefix,mean=mean,sd=sd,readlen=readlen)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)


    # step 5: report and summarize the result
    cmd= '''
    awk 'NR==FNR{{a[$1$2$3]=$4;next}}$1$2$3 in a{{print $0  a[$1$2$3]}}' {prefix}.splitsite.bed {prefix}.SR.bed >{prefix}.support.SR.bed
    cat {prefix}.support.SR.bed|awk '{{a[$1"\\t"$2"\\t"$3]=a[$1"\\t"$2"\\t"$3]","$4}}END{{for(i in a){{print i"\\t"a[i]}}}}'|sed 's/\\t,/\\t/g' > {prefix}.tmp.SR.bed
    cat {prefix}.support.PE.bed|awk '{{a[$1"\\t"$2"\\t"$3]=a[$1"\\t"$2"\\t"$3]","$4}}END{{for(i in a){{print i"\\t"a[i]}}}}'|sed 's/\\t,/\\t/g' > {prefix}.tmp.PE.bed
    awk -F $'\\t' 'BEGIN{{print "Scaffolds\\tSplit Site 1\\tSplit Site 2\\tSupporting Split Reads\\tSupporting Read Pairs(first in pair)"}}NR==FNR{{a[$1$2$3]=$4;next}}{{print $0"\\t" a[$1$2$3]}}' {prefix}.tmp.PE.bed {prefix}.tmp.SR.bed >{prefix}.tab '''.format(prefix=prefix)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)
    
    # step 6: remove intermediate file 
    cmd = '''mv {prefix}.tab {output} ; rm  {prefix}.* '''.format(prefix=prefix,output=output)
    print("****** NOW RUNNING COMMAND ******: " + cmd)
    print run_cmd(cmd)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Identify putative insertions in reference genomes from bam alignments')
    parser.add_argument('-b', '--bam', help='input bam file', required=True,dest='input',metavar='')
    parser.add_argument('-o', '--out', help='output file', required=True,dest='output',metavar='')
    parser.add_argument('-t', '--tmp', help='tmp dir', dest='tmpdir',default="tmp/",metavar='')
    parser.add_argument('-p', '--threads', help='number of threads used',type=int,default=1, required=False,dest='core',metavar='')
    parser.add_argument('-r', '--readlen', help='average read length',type=int, required=True,dest='readlen',metavar='')
    parser.add_argument('-m', '--mean', help='mean library insertion size',type=float,required=True,dest='mean',metavar='')
    parser.add_argument('-s', '--sd', help='standard deviation of library insertion size',type=float,required=True,dest='sd',metavar='')
    parser.add_argument('-n', '--number', help='number of split reads required to validate split sites',type=int, required=False,default = 4,dest='min_SR',metavar='')

    args = parser.parse_args()

    #check_tools(["sam2bed","bedtools","samtools"])
    for i in ["sam2bed","bedtools","samtools"]:
        if not is_tool(i):
            print "tool {i} is not installed".format(i=i) 
            sys.exit(0)

    pipeline(args)
