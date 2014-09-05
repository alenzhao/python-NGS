
#!/usr/bin/python
import os, glob,sys;

### program to convert indels in VCF file to full-length indels that contain the entire haplotype for the reference and indel ### 
### author Vikas Bansal, april 3 2012

def read_fasta(seqfile,fastaindex):
        sequences = [];  File = open(seqfile,'r'); lines =0; strlist = []; fastaseqs =0;
        for line in File:
                if line[0] == '>':
                        refname = line.strip().split()[0].lstrip('>');
                        sequences.append([refname,'',0]); fastaindex[refname] = fastaseqs;
                        if lines > 0:
                                #print >>sys.stderr, 'sequences',refname,fastaseqs;
                                print >>sys.stderr, 'length of list',len(strlist),sequences[fastaseqs-1][0];
                                sequences[fastaseqs-1][1] = ''.join(strlist);
                                for i in range(len(strlist)): sequences[fastaseqs-1][2] += len(strlist[i]);
                                for i in range(len(strlist)): strlist.pop();
                        fastaseqs +=1;

                else: strlist.append(line.strip().upper()); lines +=1;
                        #       print >>sys.stderr, 'length of list',len(strlist);
        sequences[fastaseqs-1][1] = ''.join(strlist);
        for i in range(len(strlist)): sequences[fastaseqs-1][2] += len(strlist[i]);
        for i in range(len(strlist)): strlist.pop();
        File.close();
        for a in range(min(fastaseqs,10)): print >>sys.stderr, sequences[a][0],sequences[a][2],sequences[a][1][0:50];
        return sequences;


def calculate_rightshift(chrom,position,ref,alt,sequences,fastaindex):

        if len(ref) > len(alt): # deletion 
                index = fastaindex[chrom];
                lengthvar = len(ref)-len(alt);
                leftpos = position; rightpos = position+lengthvar; rightshift = 0;
                print leftpos,rightpos,index,sequences[index][0],sequences[index][2];
                while sequences[index][1][leftpos] == sequences[index][1][rightpos]:  leftpos +=1; rightpos +=1; rightshift +=1;

                newref = sequences[index][1][position-1:position+lengthvar+rightshift];
                newalt = sequences[index][1][position-1]+sequences[index][1][position+lengthvar:position+lengthvar+rightshift];

                flankingL = sequences[index][1][position-20:position-1]; flankingR = sequences[index][1][position+lengthvar+rightshift:position+rightshift+20+lengthvar];

                return [newref,newalt,rightshift,flankingL,flankingR];

                #print indel[3],indel[4],
                #print '%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(chrom,position,indel[2],newref,newalt,indel[5],indel[6],indel[7],indel[8],indel[9])
        elif len(ref) < len(alt):
                index = fastaindex[chrom];
                lengthvar = len(alt)-len(ref);
                if position + 200 < sequences[index][2]: newseq = alt[1:] + sequences[index][1][position:position+200];
                else: newseq = alt[1:] + sequences[index][1][position:sequences[index][2]];
                leftpos = 0; rightpos = lengthvar; rightshift =0; maxl = len(newseq);

                while rightpos < maxl and newseq[leftpos] == newseq[rightpos]: leftpos +=1; rightpos +=1; rightshift +=1;

                newref = sequences[index][1][position-1:position+rightshift];
                newalt = alt + sequences[index][1][position:position+rightshift];
                flankingL = sequences[index][1][position-20:position-1]; flankingR = sequences[index][1][position+rightshift:position+rightshift+20];

                #print chrom,position,type,bases, sequences[index][1][position-1-rightshift:position-1],bases,rightshift;
        #       print line, chrom,position, rightshift;
                return [newref,newalt,rightshift,flankingL,flankingR];
        else:  # it is a SNP 
                index = fastaindex[chrom];
                flankingL = sequences[index][1][position-20:position-1]; flankingR = sequences[index][1][position+1:position+20];
                return [ref,alt,0,flankingL,flankingR];

def read_indels(indelfile,seqfile):
        fastaindex = {}; sequences = read_fasta(seqfile,fastaindex);

        if indelfile == '-' or indelfile == 'stdin' or indelfile == 'sys.stdin': File = sys.stdin;
        else: File = open(indelfile,'r');
        for line in File:
                if line[0] == '#':
                        #print line.strip(); 
                        continue;

                indel = line.strip('\n').split('\t'); chrom = indel[0]; position = int(indel[1]);

                ## print line for SNPs
                #if len(indel[3]) == len(indel[4]): print line.strip(); continue; 
                #if len(altalleles) ==1: continue;
                ref = indel[3]; altalleles = indel[4].split(',');

                #print line,
                for i in xrange(len(altalleles)):
                        s1=len(ref)-1; s2 = len(altalleles[i])-1;
                        if s1> s2 and s2 > 0: # both alleles are greater than 1 in length:
                                while s2 > 0 and ref[s1] == altalleles[i][s2]: s1 -=1; s2 -=1;
                                ref0 = ref[0:s1+1]; alt0 = altalleles[i][0:s2+1];
                        elif s2> s1 and s1 > 0: # both alleles are greater than 1 in length:
                                while s1 > 0 and ref[s1] == altalleles[i][s2]: s1 -=1; s2 -=1;
                                ref0 = ref[0:s1+1]; alt0 = altalleles[i][0:s2+1];
                        else: ref0 = ref; alt0 = altalleles[i];

                        #if s1 !=s2: continue;
                        if s1 == s2: continue;

                        [newref,newalt,shift,flankingL,flankingR] = calculate_rightshift(chrom,position,ref0,alt0,sequences,fastaindex);
                        #print indel[3],indel[4],ref0,alt0,shift,'|',
                        sys.stdout.write('%s\t%d\t%s\t%s\t%s\t%s\t%s\t' %(chrom,position,indel[2],newref,newalt,indel[5],indel[6]));
                        sys.stdout.write('FL=%s;OA=%s,%s;%s\t%s\t%s\n' %(flankingL + ':' + newref + ':' + flankingR,indel[3],indel[4],indel[7],indel[8],indel[9]));
                        #print flankingL + ':' + newref + ':' + flankingR;

                        #continue;

                ## take care of case where there are multiple alternate alleles...
        if indelfile == '-' or indelfile == 'stdin' or indelfile == 'sys.stdin': pass;
        else: File.close();



if len(sys.argv) < 3: print  'python leftjustify-indels.py indel-file reference.fasta'; sys.exit();
read_indels(sys.argv[1],sys.argv[2]);



