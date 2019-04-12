import pandas as pd
import numpy as np
import os
from tqdm import tqdm

if not os.path.isfile('data/GencodeGenes_v19_hg19.gz'):
    genes = pd.read_csv('data/gencode_v19_hg19.gz', compression='gzip', sep="\t")
    transcript_mapping_file = pd.read_csv('data/ENST_ENSG_mapping_file.txt', sep='\t')
    transcript_names = [gene.split('.')[0] for gene in genes['#name']]
    rowIndices = [ list(transcript_mapping_file['Transcript stable ID']).index(transcript_name) for transcript_name in tqdm(transcript_names) ]
    ENSG_list = list(transcript_mapping_file.iloc[ rowIndices,: ]['Gene stable ID'])
    genes['name3'] = ENSG_list
    genes.to_csv('data/GencodeGenes_v19_hg19.gz', compression='gzip', sep="\t", encoding='utf-8', index=False)
else:
    genes = pd.read_csv('data/GencodeGenes_v19_hg19.gz', compression='gzip', sep="\t", encoding='utf-8')


# Can't get collapse_annotation.py to work, so here's my short version for collapsing the genes:
class Exon:
    def __init__(self, number, start_pos, end_pos):
        self.number = int(number)
        self.start_pos = start_pos
        self.end_pos = end_pos

class Transcript:
    def __init__(self, transcript_num, start_pos, end_pos, exons):
        self.transcript_num = int(transcript_num)
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.exons = exons

class Gene:
    def __init__(self, ensg_name, gene_name, chrom, strand, start_pos, end_pos, transcripts):
        self.ensg_name = ensg_name
        self.name = gene_name
        self.chrom = chrom
        self.strand = strand
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.transcripts = transcripts

def initializeBins(transcript):
    '''
    Returns start and end positions of exons in a dataframe
    '''
    starts = []
    ends = []
    for exon in transcript.exons:
        starts.append(exon.start_pos)
        ends.append(exon.end_pos)
    bins = pd.DataFrame({'starts': starts, 'ends': ends})
    return bins

def cleanup(bins):
    sorted_bins = bins.sort_values(by=['starts'])
    cleanedBins = pd.DataFrame(columns=['starts','ends'])
    i=0
    while i < (sorted_bins.shape[0]-1):
        currBin = list(sorted_bins.iloc[i,:])
        nextBin = list(sorted_bins.iloc[i+1,:])
        while overlap(currBin, nextBin) and i < (sorted_bins.shape[0]-1):
            currBin = overlapExtend(currBin, nextBin)
            i += 1
            try:
                nextBin = list(sorted_bins.iloc[i+1,:])
            except:
                pass
        cleanedBins = cleanedBins.append(pd.DataFrame({'starts':[currBin[0]], 'ends':[currBin[1]]}), ignore_index=True)
        i += 1
    if i != sorted_bins.shape[0]:
        i += 1
        cleanedBins = cleanedBins.append(pd.DataFrame({'starts':[nextBin[0]], 'ends':[nextBin[1]]}), ignore_index=True)
    return cleanedBins

def createTranscript(bins):
    bins = cleanup(bins)
    starts = list(bins['starts'])
    ends = list(bins['ends'])
    minStart = starts[0]
    maxEnd = ends[-1]
    exons = []
    for i in np.arange(len(starts)):
        exons.append(Exon(i, starts[i], ends[i]))
    return Transcript(i, minStart, maxEnd, exons)

def overlap(bin1, bin2):
    '''
    bin1 and bin2 are lists of length 2 each
    Returns True if there is overlap, False otherwise
    '''
    return (bin2[0]>=bin1[0] and bin2[0]<=bin1[1]) or (bin2[1]>=bin1[0] and bin2[1]<=bin1[1])

def overlapExtend(bin1, bin2):
    '''
    bin1 and bin2 must overlap and are lists of length 2 each
    Returns a new bin with min and max of the extended bin
    '''
    allpositions = []
    allpositions.extend(bin1)
    allpositions.extend(bin2)
    minBin = min(allpositions)
    maxBin = max(allpositions)
    newBin = [minBin, maxBin]
    return newBin

def collapse(gene):
    '''
    Modifies gene.transcripts to be one Transcript with "unionized" exons
    '''
    if len(gene.transcripts) == 1:
        return gene
    else:
        bins = initializeBins(gene.transcripts[0])
        for transcript in gene.transcripts[1:]:
            for exon in transcript.exons:
                overlappingExon = False
                exonBin = [exon.start_pos, exon.end_pos]
                for i in np.arange(bins.shape[0]):
                    currBin = list(bins.iloc[i,:])
                    if overlap(currBin, exonBin):
                        currBin = overlapExtend(currBin, exonBin)
                        bins.iloc[i,:]['starts'] = currBin[0]
                        bins.iloc[i,:]['ends'] = currBin[1]
                        overlappingExon = True
                if not overlappingExon:
                    bins = bins.append(pd.DataFrame({'starts':[exon.start_pos], 'ends':[exon.end_pos]}), ignore_index=True)
        gene.transcripts = [createTranscript(bins)]
        return gene

gene_groups = genes.groupby('name3')
all_genes_collapsed = []
for gene in tqdm(gene_groups): # for each gene
    ensg_name = gene[0]
    gene_df = gene[1]
    gene_name = list(gene_df['name2'])[0]
    chrom = list(gene_df['chrom'])[0]
    strand = list(gene_df['strand'])[0]
    transcripts = []
    minStart = 0
    maxEnd = 0
    for i in np.arange(gene_df.shape[0]): # for each transcript
        starts = [int(num.strip()) for num in list(gene_df['exonStarts'])[i].split(',') if num.strip() != ""]
        ends = [int(num.strip()) for num in list(gene_df['exonEnds'])[i].split(',') if num.strip() != ""]
        if minStart==0:
            minStart = min(starts)
        if maxEnd==0:
            maxEnd = max(ends)
        if min(starts) < minStart:
            minStart = min(starts)
        if max(ends) > maxEnd:
            maxEnd = max(ends)
        exons = []
        for j in np.arange(len(starts)): # for each exon
            exons.append(Exon(j, starts[j], ends[j]))
        transcripts.append(Transcript(i, min(starts), max(ends), exons))
    currGene = Gene(ensg_name, gene_name, chrom, strand, minStart, maxEnd, transcripts)
    all_genes_collapsed.append(collapse(currGene))

