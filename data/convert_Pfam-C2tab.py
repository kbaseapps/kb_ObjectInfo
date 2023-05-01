
path = '/Users/mland/my_KBASE/kb_ObjectInfo/data/'

cat = 'Other'

outfile = open(path + "Pfam-C.clans.v35.tsv","w")

with open(path + "Pfam-C.tsv", "r") as file1:
    for line in file1:
        newline = line.replace("\n", "")
        (header,linetype,desc) = newline.split("\t")[0:3]
        
        if header == '# ST' and cat != 'Other':
            outline = "\t".join([cat,cat_id,cat_name,"\n"])
            outfile.write(outline)
            cat = 'Other'
            cat_name = 'Other'
            cat_group = 'Other'
            cat_id    = 'Other'
        elif linetype == "AC":
            cat = desc
            if '.' in cat and len(cat) < 15:
                cat = cat.split('.')[0]
        elif linetype == "ID":
            cat_id = desc
        elif linetype == "DE":
            cat_name = desc
            cat_name = cat_name.replace('"','')

file1.close()
outfile.close()
