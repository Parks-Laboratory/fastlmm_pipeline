# -*- coding: utf-8 -*-
"""
Retrieves Mouse Diversity Array genotypes from database in PLINK format
"""

import pyodbc
import os
import sys

def get_genotypes(strains, iids, output_fn):
    query_template = '''select snp_chr, rsID, 0 as centimorgans, snp_bp_mm10, %s
    from HMDP.dbo.genotype_calls_plink_format
    order by snp_chr, snp_bp_mm10
    '''

    outfile = open(output_fn, 'w')
    c = pyodbc.connect('SERVER=PARKSLAB;DATABASE=HMDP;Trusted_Connection=Yes;DRIVER={SQL Server}')
    
    ori = c.execute('select top 1 * from HMDP.dbo.genotype_calls_plink_format')
    for row in ori:
        cols = [t[0] for t in row.cursor_description]

    q = query_template % ', '.join(['[%s]' % x for x in strains if x in cols])
    res = c.execute(q)
    tfam = 0
    linebuffer = []
    # maybe retrieve all rows, then print at end
    for row in res:
        if tfam == 0:
            # generate .tfam file for PLINK
            colnames = [t[0] for t in row.cursor_description]
            # sanitize strain names
            colnames = [x.replace('/', '.').replace(' ', '.') for x in colnames]
            #if len(set(colnames)) < len(colnames):
            #    # need to make identifiers unique
            tfam = 1
            tfam_outfile = open(output_fn.replace('.tped', '.tfam'), 'w')
            for i, fid in enumerate(colnames[4:]):
                iid = iids[i].replace('/', '.').replace(' ', '.')
                tfam_outfile.write( '\t'.join(map(str, [fid, iid, 0, 0, 0, -9])) )
                tfam_outfile.write('\n')
            tfam_outfile.close()


        # linebuffer.append('\t'.join(map(str, row)))
        # convert missing data, join line to format
        linebuffer.append('\t'.join( map(str, ['0 0' if x is None else x for x in list(row)]) ))
        # print 50000 lines at a time
        if len(linebuffer) == 50000:
            outfile.write('\n')
            outfile.write( '\n'.join(linebuffer) )
            outfile.flush()
            linebuffer = []

    # flush lines if any are left over
    if linebuffer:
        outfile.write('\n')
        outfile.write( '\n'.join(linebuffer) )

    c.close()
    outfile.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        lines = [x.strip().split('\t') for x in sys.stdin.readlines()]
        
        strains, iids = zip(*lines)
        print(strains)
        get_genotypes(strains, iids, 'test.tped')
    else:
        for filename in sys.argv[1:]:
            print('Processing' + filename)
            f = open(filename)
            lines = [x.strip().split('\t') for x in f.readlines() if x.strip()]
            f.close()

            strains, iids = zip(*lines)
            output_fn = '%s.tped' % os.path.splitext(filename)[0]
            get_genotypes(strains, iids, output_fn)
            print('Genotypes saved to' + output_fn)
