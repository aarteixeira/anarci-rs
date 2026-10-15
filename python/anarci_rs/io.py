"""ANARCI text, CSV and FASTA compatibility helpers.

Adapted from ANARCI, copyright (C) 2016 Oxford Protein Informatics Group.
Distributed under the BSD 3-Clause License; see licenses/ANARCI-BSD-3-Clause.txt.
"""
import gzip
from itertools import groupby, islice
from textwrap import wrap

amino_acids = sorted("QWERTYIPASDFGHKLCVNM")
set_amino_acids = set(amino_acids)
chain_type_to_class = {"H":"H", "K":"L", "L":"L", "A":"A", "B":"B", "G":"G", "D":"D"}

def read_fasta(filename):
    """Read plain or gzip FASTA records, preserving full header text."""
    return list(fasta_iter(filename))


def fasta_iter(fasta_name):
    """Yield (header, sequence) pairs from a plain or gzip FASTA file."""
    filename = str(fasta_name)
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(filename, "rt") as handle:
        groups = (group for _, group in groupby(handle, lambda line: line.startswith(">")))
        for header in groups:
            name = next(header)[1:].strip()
            sequence = "".join(line.strip() for line in next(groups, ()))
            yield name, sequence


def write_fasta(sequences, f):
    """
    Write a list of sequences to file. 
    
    should be a list of name, sequence tuples
    
    f should be an open file
    """
    for name, sequence in sequences:
        print(">%s"%name, file=f)
        print('\n'.join(['\n'.join(wrap(block, width=80)) for block in sequence.splitlines()]), file=f)

def validate_sequence(sequence):
    """
    Validate the length and canonical amino acid alphabet of a protein sequence.
    """
    assert len(sequence) < 10000, "Sequence too long."
    assert not (set( sequence.upper() ) - set_amino_acids ), "Unknown amino acid letter found in sequence: %s"% ", ".join(list((set( sequence.upper() ) - set_amino_acids )))    
    return True

def validate_numbering(xxx_todo_changeme, name_seq=[]):
    """
    Wrapper to do some basic validation of the numbering.
    
    Further validation could be done but at the moment we just check that the numbering indices are incremental (they should be)
    """
    (numbering, start, end) = xxx_todo_changeme
    name, seq = name_seq
    last = -1
    nseq=""

    for (index, _), a in numbering:
        assert index >= last, "Numbering was found to decrease along the sequence %s. Please report."%name
        last = index
        nseq += a.replace("-","")

    assert nseq in seq.replace("-",""), "The algorithm did not number a contiguous segment for sequence %s. Please report"%name

    return numbering, start, end

def grouper(n, iterable):
    '''
    Group entries of an iterable by n
    '''
    it = iter(iterable)
    def take():
        while 1:
            yield list( islice(it,n) )
    return iter(take().__next__, [] )

def anarci_output(numbered, sequences, alignment_details, outfile, sequence_id=None, domain_id=None):
    """
    Outputs to open file

    If sequence_id is specified as an integer then only this sequence will be printed. 
    Otherwise all sequences will be printed.

    If domain_id is specified as an integer then only this domain will be printed. 
    Otherwise all domains will be printed.

    If domain_id is specified then sequence_id must also be specified. 
    """       
    assert (sequence_id is not None) or (sequence_id is None and domain_id is None), "If domain_id is specified, sequence_id must also be specified."
    for i in range(len(numbered)):
        if sequence_id is None:
            print("# %s"%sequences[i][0], file=outfile) # print the name
        if numbered[i] is not None:
            if sequence_id is not None:
                if i != sequence_id: continue
            print("# ANARCI numbered", file=outfile)
            for j in range( len(numbered[i])): # Iterate over domains
                if domain_id is not None:
                    if j != domain_id: continue
                print("# Domain %d of %d"%(j+1, len(numbered[i]) ), file=outfile)
                print("# Most significant HMM hit", file=outfile)
                print("#|species|chain_type|e-value|score|seqstart_index|seqend_index|", file=outfile)
                alignment_details[i][j]["evalue"] = str( alignment_details[i][j]["evalue"] )
                print("#|%s|%s|%s|%.1f|%d|%d|"%tuple( [alignment_details[i][j][field] for field in 
                                                                     ["species","chain_type","evalue","bitscore"]] 
                                                                   +[ numbered[i][j][1], numbered[i][j][2]] ), file=outfile)
                
                if 'germlines' in alignment_details[i][j]:
                    print('# Most sequence-identical germlines', file=outfile)
                    print('#|species|v_gene|v_identity|j_gene|j_identity|', file=outfile)
                    (species, vgene), vid =alignment_details[i][j]['germlines'].get('v_gene', [['','unknown'],0])
                    if vgene is None:
                        vgene, vid = 'unknown', 0
                    (_,jgene), jid =alignment_details[i][j]['germlines'].get('j_gene', [['','unknown'],0])
                    if jgene is None:
                        jgene, jid = 'unknown', 0
                    print('#|%s|%s|%.2f|%s|%.2f|'%(species, vgene, vid, jgene, jid ), file=outfile)	
                chain_type = chain_type_to_class[  alignment_details[i][j]["chain_type"] ]
                print("# Scheme = %s"%alignment_details[i][j]["scheme"], file=outfile)
                if len( numbered[i][j][0] ) == 0:
                    print("# Warning: %s scheme could not be applied to this sequence."%alignment_details[i][j]["scheme"], file=outfile)            
                for (index, insertion), aa in numbered[i][j][0]:
                    print(chain_type, ("%d"%index).ljust(5), insertion, aa, file=outfile)
        print("//", file=outfile)

def csv_output(sequences, numbered, details, outfileroot):
    '''
    Write numbered sequences to csv files. A csv file is written for each chain type.

    Kappa and Lambda chains are written to the same file

    The sequences will written aligned to the numbering scheme. Gaps in the sequences with respect to the alignment are written
    as a '-'

    @param sequences: List of name, sequence tuples    
    @param numbered: Numbered sequences in the same order as the sequences list. 
    @param details: List of alignment details in the same order as the sequences list.
    @param outfileroot: The file path for csv files to write. _<chain_type>.csv will be appended to this.
    '''
    
    chain_types = {}        
    pos_ranks = {}
    all_pos = {}
    _lc = {'K':'KL','L':'KL'}


    # Divide the set into chain types and find how to order the numbering for each type.
    for i in range( len(sequences) ): # Iterate over entries
        if numbered[i] is None: continue

        for j in range(len(numbered[i])): # Iterate over domains.
            # Record the chain type index
            c = details[i][j]['chain_type']
            c = _lc.get(c, c) # Consider lambda and kappa together.
            chain_types.setdefault( c, [] ).append( (i,j) ) 
            if c not in pos_ranks:
                pos_ranks[c] = {}
                all_pos[c] = set()

            # Update the insertion order for the scheme. i.e. is it A B C or C B A (e.g. imgt 111 and 112 repectively)
            l = -1 
            r = 0 
            for p, _ in numbered[i][j][0]:
                if p[0] != l:
                    l = p[0]
                    r = 0
                else:
                    r +=1
                pos_ranks[c][p] = max( r, pos_ranks[c].get( p, r ) )
                all_pos[c].add( p )

    # Write a new file for each chain type. Kappa and lambda are written together as light chains. 
    for cts in ['H','KL','A','B','G','D']:
        if cts in chain_types:
            with open( outfileroot + '_%s.csv'%cts, 'w' ) as out:

                # Sort the positions by index and insertion order
                positions = sorted( all_pos[cts], key = lambda p: (p[0], pos_ranks[cts][p]) )

                # Header line
                fields = ['Id','domain_no','hmm_species','chain_type','e-value','score','seqstart_index','seqend_index',
                          'identity_species','v_gene','v_identity','j_gene','j_identity']
                fields += [ ('%d%s'%(p)).strip() for p in positions ]
                print(','.join( fields ), file=out)

                # Iterate over the domains identified
                for i,j in chain_types[cts]:
                    line = [ sequences[i][0].replace(',',' '),      
                             str(j),
                             details[i][j].get('species',''),
                             details[i][j].get('chain_type',''),
                             str(details[i][j].get('evalue','')),
                             str(details[i][j].get('bitscore','')),
                             str(numbered[i][j][1]),
                             str(numbered[i][j][2]),
                             details[i][j].get('germlines',{}).get( 'v_gene',[['',''],0] )[0][0],
                             details[i][j].get('germlines',{}).get( 'v_gene',[['',''],0] )[0][1],
                             '%.2f'%details[i][j].get('germlines',{}).get( 'v_gene',[['',''],0] )[1],
                             details[i][j].get('germlines',{}).get( 'j_gene',[['',''],0] )[0][1],
                             '%.2f'%details[i][j].get('germlines',{}).get( 'j_gene',[['',''],0] )[1] ]

                    # Hash the numbering. Insertion order has been preserved in the positions sort.
                    d = dict( numbered[i][j][0] )
                    line += [ d.get(p,'-') for p in positions ]

                    assert len( line ) == len( fields )
                    print(','.join( line ), file=out)
