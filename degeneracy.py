#!/usr/local/bin/python3
import sys
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

matchParent = re.compile('Parent=([^;]+)')
uniqCheck = set()

def define_sites():
	Degenerate = {
		"4fold":{"TCT","TCC","TCG","TCA","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
			"CGT","CGC","CGA","CGG","ACT","ACC","ACA","ACG","GTT","GTC","GTA","GTG",
			"GCT","GCC","GCA","GCG","GGT","GGC","GGA","GGG"
		},
		
		"2fold":{"CAG", "CAA", "GAG", "GAA", "AAG", "AAA", "TTG", "TTA", "AGG", "AGA", "CAT", "CAC", "GAT", "GAC", "AAT", "AAC", "TAT", "TAC", "TTT", "TTC", "AGT", "AGC", "TGT", "TGC"
		},

        "0fold_1":{"TGA","TAG","TAA"
		},
		
		"0fold_2":{"AGG","AGA","CGA","CGG","CGC","CGT","CTG","CTC","CTA","CTT","TTG","TTA"
		},

		"0fold_12":{"GGA","GGC","GGG","GGT","GAG","GAA","GAT","GAC","GCA","GCC","GCA",
			"GCT","GTA","GTG","GTC","GTT","AAG","AAA","AAC","AAT","ACG","ACC","ACT",
			"ACA","ATA","ATC","ATT","CAG","CAA","CAT","CAC","CCG","CCC","CCA","CCT",
			"TGC","TGT","TAC","TAT","TTC","TTT"
		},

		"0fold_all":{"ATG","TGG"
		}

	}

	return(Degenerate)


def main():
    #verify input arguments
    argumentCounts =  len(sys.argv)
    if (argumentCounts!=3):
        print("Error: Two arguments are required!");
        help_menu()
        sys.exit()

    gffFileName = sys.argv[1]
    genomeFilename =  sys.argv[2]

    if (not os.path.isfile(gffFileName)):
        print(f"Error: File {gffFileName} does not exist!")
        help_menu()
        sys.exit()
    if (not os.path.isfile(genomeFilename)):
        print(f"Error: File {genomeFilename} does not exist!")
        help_menu()
        sys.exit()

    chrToSeq = SeqIO.to_dict(SeqIO.parse(genomeFilename, "fasta"))


    #parse the gff file
    transcriptToChrStrand={}
    transcriptToCDS={}
    with open(gffFileName, 'r') as gffh:
        for line in gffh:
            if (line.startswith("#")):
                continue
            if (not re.search('\d', line)):
                continue
            fieldArray = line.split(sep="\t")
            chrName = fieldArray[0]
            featureSource = fieldArray[1]
            featureType = fieldArray[2]
            startPos = int(fieldArray[3])
            endPos = int(fieldArray[4])
            strand = fieldArray[6]
            gffphase= fieldArray[7]
            geneAnnot = fieldArray[8]

            #geneId = getParent(geneAnnot)

            #print(geneId)
            if (featureType=="CDS"):
                parentTranscript = getParent(geneAnnot)
                if parentTranscript in transcriptToCDS:
                    transcriptToCDS[parentTranscript].append((startPos,endPos))
                else:
                    transcriptToChrStrand[parentTranscript] =(chrName,strand)
                    transcriptToCDS[parentTranscript] = [(startPos,endPos)]                 
            else:
                pass
    gffh.close()

    # process each transcript, record CDS range, output CDS file and protein file
    degenerate = define_sites()
    myCDSput = open("cds.fasta", "w")
    myPTput = open("protein.fasta", "w")
    myDEG0 = open("deglist_0fold.bed", "w")
    myDEG4 = open("deglist_4fold.bed", "w")
    myDEG2 = open("deglist_2fold.bed", "w")


    for mRNAID in transcriptToCDS:

        chrName = transcriptToChrStrand[mRNAID][0]
        if (not chrName in chrToSeq):
            print ("Error: no chr named " + chrName)
            continue;
        strand = transcriptToChrStrand[mRNAID][1]
        reverseflag=False
        if strand=="-":
            reverseflag = True
        transcriptToCDS[mRNAID].sort(key=lambda x: x[0], reverse=reverseflag)

        CDSSeq = Seq("", generic_dna)
        descriptionLine = chrName + "(" + strand + ") Segments:"

        CDSPosToChrPos = []
        for rangecoord in transcriptToCDS[mRNAID]:
            cds_start = rangecoord[0]-1
            cds_end = rangecoord[1]
            descriptionLine+=str(rangecoord[0]) + "-" + str(rangecoord[1]) + ";"
            cds_segment = chrToSeq[chrName].seq[cds_start:cds_end];

            if reverseflag:  ## minus strand
                cds_segment = cds_segment.reverse_complement()
                CDSPosToChrPos.extend(range(cds_end -1,cds_start-1,-1))
            else:   ## plus strand
                CDSPosToChrPos.extend(range(cds_start,cds_end))

            CDSSeq += cds_segment
        
        # output the CDS and translation
        descriptionLine = descriptionLine.rstrip(';')
        SeqIO.write(SeqRecord(CDSSeq, id=mRNAID, description=descriptionLine), myCDSput, "fasta")
        SeqIO.write(SeqRecord(CDSSeq.translate(), id=mRNAID, description=descriptionLine), myPTput, "fasta")

        #identify all degenerate codons and put in a list
        degList = []
        for cdsPosIndex in range(0, len(CDSSeq), 3):
            tri_nt = str(CDSSeq[cdsPosIndex : (cdsPosIndex+3)])


            # 0 fold
            typestr = "" 
            codonPosIndex = CDSPosToChrPos[cdsPosIndex]
            wobblePosOffsets = []
            if tri_nt in degenerate["0fold_1"]:
                typestr="0fold_1"
                wobblePosOffsets=[0]
            elif tri_nt in degenerate["0fold_2"]:
                typestr="0fold_2"
                wobblePosOffsets=[1]
            elif tri_nt in degenerate["0fold_12"]:
                typestr="0fold_12"
                wobblePosOffsets=[0,1]
            elif tri_nt in degenerate["0fold_all"]:
                typestr="0fold_all"
                wobblePosOffsets=[0,1,2]
            
            if typestr != "":
                for offsetValue in wobblePosOffsets:
                    wobblePosIndex = CDSPosToChrPos[cdsPosIndex +offsetValue]
                    write_site(myDEG0, typestr, chrName, strand, mRNAID, codonPosIndex, cdsPosIndex, offsetValue, wobblePosIndex, tri_nt)

            ## 2fold and 4fold
            typestr = "" 
            if tri_nt in degenerate["4fold"]:
                typestr="4fold"
            elif tri_nt in degenerate["2fold"]:
                typestr="2fold"
            else:
                continue

            codonPosIndex = CDSPosToChrPos[cdsPosIndex]
            wobblePosIndex = CDSPosToChrPos[cdsPosIndex +2]
            offsetValue = 2

            if typestr=="4fold":
                FH = myDEG4
            elif typestr=="2fold":
                FH = myDEG2

            write_site(FH, typestr, chrName, strand, mRNAID, codonPosIndex, cdsPosIndex, offsetValue, wobblePosIndex, tri_nt)

#            myID= chrName + ":" + str(wobblePosIndex)
#            if myID in uniqCheck:
#                continue
#            else:
#                uniqCheck.add(myID)
#            namestr = ":".join([typestr, mRNAID, str(cdsPosIndex), str(codonPosIndex), tri_nt])
#
#            if typestr=="4fold":
#                myDEG4.write("\t".join([chrName, str(wobblePosIndex), str(wobblePosIndex+1), namestr, ".", strand]))
#                myDEG4.write("\n")
#            elif typestr=="2fold":
#                myDEG2.write("\t".join([chrName, str(wobblePosIndex), str(wobblePosIndex+1), namestr, ".", strand]))
#                myDEG2.write("\n")

# verification code
#                if reverseflag:
#                    chrnt = chrToSeq[chrName].seq[chrPosIndex:chrPosIndex+3].reverse_complement();
#                else:
#                    chrnt = chrToSeq[chrName].seq[chrPosIndex-2:chrPosIndex+1];
#                
#                degList.append((cdsPosIndex, chrPosIndex, tri_nt, str(chrnt)))
                
#                 if (tri_nt !=str(chrnt)):
#                    #print(cdsPosIndex, chrPosIndex, CDSPosToChrPos[cdsPosIndex], tri_nt, str(chrnt), descriptionLine)
#                    pass


    myCDSput.close()
    myPTput.close()
    myDEG2.close()
    myDEG4.close()
    myDEG0.close()

    os.system("sort -k1,1V -k2,2n deglist_2fold.bed > deglist_2fold_sorted.bed")
    os.system("sort -k1,1V -k2,2n deglist_4fold.bed > deglist_4fold_sorted.bed")
    os.system("sort -k1,1V -k2,2n deglist_0fold.bed > deglist_0fold_sorted.bed")
    #sort -k1,1V -k2,2n deglist_2fold.bed > deglist_2fold_sorted.bed


##
def write_site(FileHandle, typestr, chrName, strand, mRNAID, codonPosIndex, cdsPosIndex, offsetValue, wobblePosIndex, tri_nt):
    myID= chrName + ":" + str(wobblePosIndex)
    if myID in uniqCheck:
        return
    else:
        uniqCheck.add(myID)
    namestr = ":".join([typestr, mRNAID, str(cdsPosIndex), str(codonPosIndex), str(offsetValue), tri_nt])
    FileHandle.write("\t".join([chrName, str(wobblePosIndex), str(wobblePosIndex+1), namestr, ".", strand]))
    FileHandle.write("\n")
    return


def help_menu():
    commandStr = sys.argv[0]
    helpstr = f"""
Usage:
{commandStr} annotation.gff3 genome.fasta
    """
    print(helpstr)


def getParent(geneStr):
    geneId = matchParent.search(geneStr).group(1)
    return geneId

if __name__ == '__main__':    
    main()


#sort -k1,1V -k2,2n deglist_2fold.bed > deglist_2fold_sorted.bed
#sort -k1,1V -k2,2n deglist_2fold.bed > deglist_2fold_sorted.bed