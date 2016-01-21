class fasta:

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(fileFASTA):
        '''Checks for fasta by file extension'''
        file_lower = fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna'):
            with open(fileFASTA, "r") as f:
                return ParseFASTA(f)

    def ParseFASTA(fileFasta):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFasta:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        return fasta_list
