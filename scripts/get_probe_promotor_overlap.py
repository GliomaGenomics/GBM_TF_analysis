
probes={}
with open('downloaded_data/infinium-methylationepic-v-1-0-b5-manifest-file_extract.txt', 'r') as probesfile:
    probesfile.readline()
    probesfile.readline()
    probesfile.readline()
    probesfile.readline()
    probesfile.readline()
    probesfile.readline()
    probesfile.readline()
    probesfile.readline()
    for line in probesfile:
        l=line.strip().split()
        if l[0].startswith('cg'):
            if l[1] not in probes:
                probes[l[1]]=[]
            probes[l[1]].append([l[2],l[3],l[0],[]])

with open('gsea_files/gencode.v27.chr_patch_hapl_scaff.annotation_PromotersTSS_1000.txt','r') as promotorfile:
    promotorfile.readline()
    for line in promotorfile:
        l=line.split()
        if l[0] in probes:
            for probe in probes[l[0]]:
                if int(probe[0])+1>=int(l[1]) and int(probe[0])+1<=int(l[2]) and l[4][0:15] not in probe[3]:
                    probe[3].append(l[4][0:15])             

with open('methylation/probe_promotor_overlap.txt', 'w+') as outfile:
    for c in probes:
        for p in probes[c]:
            pp=[c]
            pp.extend(p[:3])
            pp.append(','.join(p[3]))
            outfile.write('\t'.join(pp))
       	    outfile.write('\n')

