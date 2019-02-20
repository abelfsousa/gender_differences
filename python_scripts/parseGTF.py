import sys

def readGTF(gtfFile):
    l = []
    r = open(gtfFile, "r")
    for line in r:
        if line[0] != "#":
            line = line.strip().split("\t")
            if line[2] == "gene":
                chrom = line[0]
                atts = line[8].split(";")
                del atts[-1]
                gene_id = atts[0].split(" ")[1].replace('"','')
                gene_type = atts[1].split(" ")[2].replace('"','')
                gene_name = atts[3].split(" ")[2].replace('"','')
                l.append((gene_id, gene_name, gene_type, chrom))
    r.close()
    return l


def writePositions(refGTF):
    GTFlist = readGTF(refGTF)
    w = open("geneAnnot.gencode.v22.txt", "w")
    w.write("%s\t%s\t%s\t%s\n" % ("geneID", "geneName", "geneType", "chrom"))
    for gene in GTFlist:
        w.write("%s\t%s\t%s\t%s\n" % (gene[0], gene[1], gene[2], gene[3]))
    w.close()


GTFfile = sys.argv[1]
writePositions(GTFfile)