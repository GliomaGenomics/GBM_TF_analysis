with open('sequences/le70_sequence_positions.txt', 'r') as infile:
    with open('sequences/le70_sequence_positions_merged.txt', 'w+') as outfile:
        active=infile.readline().split()
        for line in infile:
            l=line.split()
            if l[0]==active[0]:
                if int(l[2]) < int(active[3]):
                    active[3]=l[3]
                else:
                    outfile.write('\t'.join(active)+'\n')
                    active=l
            else:
                outfile.write('\t'.join(active)+'\n')
                active=l
        outfile.write('\t'.join(active)+'\n')

with open('sequences/control_sequence_positions.txt', 'r') as infile:
    with open('sequences/control_sequence_positions_merged.txt', 'w+') as outfile:
        active=infile.readline().split()
        for line in infile:
            l=line.split()
            if l[0]==active[0]:
                if int(l[2]) < int(active[3]):
                    active[3]=l[3]
                else:
                    outfile.write('\t'.join(active)+'\n')
                    active=l
            else:
                outfile.write('\t'.join(active)+'\n')
                active=l
        outfile.write('\t'.join(active)+'\n')

