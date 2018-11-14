import sqlite3
import xlsxwriter
from ngsscriptlibrary import MetricsDBReader
from ngsscriptlibrary import get_rs_gpos_dict

def write_list(listname, worksheet, row=0, col=0, skip=1, header=False,
               orientation='rows', format=None, formatheader=None):
    if header:
        worksheet.write(row, col, header, formatheader)
        row += skip

    if orientation == 'rows':
        [worksheet.write(row + i, col, ii, format)
         for i, ii in enumerate(listname)]
        row = row + len(listname) + 2

    elif orientation == 'cols':
        [worksheet.write(row, col + i, ii, format)
         for i, ii in enumerate(listname)]
        col = col + len(listname) + 2
        row += 1

    return (row, col)


class CreateReport:

    def __init__(self, sample=None, serie=None, db_metrics=None, db_samplesheet=None, db_targets=None):
        if db_samplesheet is None:
            db_samplesheet = 'samplesheets.sqlite'
        if db_metrics is None:
            db_metrics = 'metrics.sqlite'
        if db_targets is None:
            db_targets = 'captures.sqlite'            
        if sample is None and serie is None:
            raise ValueError('Geen sample of serie opgegeven')
        elif sample is not None and serie is None:
            sql = "SELECT * FROM todo WHERE SAMPLE='{}'".format(sample)
        elif sample is None and serie is not None:
            sql = "SELECT * FROM todo WHERE SERIE='{}'".format(serie)
        elif sample is not None and serie is not None:
            sql = "SELECT * FROM todo WHERE (SERIE='{}' AND SAMPLE='{}')".format(serie, sample)
        
        conn = sqlite3.connect(db_samplesheet)
        c = conn.cursor()
        c.execute(sql)
        
        self.todo = dict()
        self.db_metrics = db_metrics
        
        for _ in c.fetchall():
            serie, sample, genesis, capture, pakket, panel, cnv_screen, cnv_diag, mosa = _
            self.todo[sample] = dict()
            self.todo[sample]['serie'] = serie
            self.todo[sample]['genesis'] = genesis
            self.todo[sample]['capture'] = capture
            self.todo[sample]['pakket'] = pakket
            self.todo[sample]['panel'] = panel
            self.todo[sample]['cnv_screen'] = bool(cnv_screen)
            self.todo[sample]['cnv_diag'] = bool(cnv_diag)
            self.todo[sample]['mosa'] = bool(mosa)
            
        conn.close()
        
        conn = sqlite3.connect(db_targets)
        c = conn.cursor()

        for sample in self.todo:
            cap, v = self.todo[sample]['capture'].split('v')
            c.execute("SELECT oid FROM captures WHERE (capture='{}' AND versie='{}')".format(cap, v))
            oid = c.fetchone()[0]
            self.todo[sample]['oid'] = oid
        
        conn.close()
            
        
        self.conn = sqlite3.connect(db_metrics)
        self.c = self.conn.cursor()
            
    def write_excel(self, sample=None):
    
        if sample is None and not self.sample:
            raise ValueError('Geen sample opgegeven')
        elif sample is None:            
            sample = self.sample
    
        REF = 'lifescope.hg19.fa'
        DBSNP = 'dbsnp137.hg19.vcf'
        GATK = 'GATK3.8 HaplotypeCaller'
        PICARD = 'picard-tools 1.95'
        BWA = 'bwa-mem 0.7.12-r1039'

        serie = self.todo[sample]['serie']
        capture = self.todo[sample]['capture']            
    
        wb = xlsxwriter.Workbook('{}.xlsx'.format(sample))
        wb.set_properties({
            'title':    sample,
            'subject':  'MiSEQUENCING',
            'author':   'Martin Haagmans',
            'comments': 'Created with Python and XlsxWriter'})
            
        ws1 = wb.add_worksheet('patient + Miseq info')
        ws2 = wb.add_worksheet('low coverage')
        ws3 = wb.add_worksheet('varianten')

        headerformat = wb.add_format()
        headerformat.set_font_size(16)
        
        underlined = wb.add_format()
        underlined.set_bottom()
        
        gray = wb.add_format()
        gray.set_bg_color('gray')

        ws1.set_column('A:A', 35)

        ws2.set_column('A:A', 25)
        ws2.set_column('F:H', 13)

        ws3.set_column('A:C', 16)
        ws3.set_column('D:D', 8)
        ws3.set_column('E:E', 16)
        ws3.set_column('F:F', 11)
        ws3.set_column('G:G', 13)
        ws3.set_column('H:H', 22)
        ws3.set_column('I:I', 14)
        ws3.set_column('J:J', 18)
        ws3.set_column('K:K', 13)
        ws3.set_column('L:L', 25)

        MISEQ = ['NGS DNA nummer:', 'Sanger DNA nummer:', 'Familienummer:',
                 'Geboortedatum:', 'Geslacht:', 'Diagnose:', 'Samenvatting:']

        SEQRAP = ['Reads:', '% chimeric reads:',
                  'Gemiddelde lengte gemapte forward:',
                  '% gemapte forward paired:',
                  'Gemiddelde lengte gemapte reverse:',
                  '% gemapte reverse paired:', 'Target (bp):',
                  '% unieke reads:', '% ontarget:', '',
                  'Informatieve coverage:', 'Standaarddeviatie:']

        SNPCHECK = ['Locus', 'rsID', 'NGS', 'TaqMan', 'Result', 'Paraaf staf']

        INFO = ['Sample', 'Pakket', 'Panel', 'OID capture', 'Picard', 'GATK',
                'Recalibrate dbSNP', 'Referentie', 'Aligner']

        INFODATA = ['{}'.format(sample),
                    '{}'.format(self.todo[sample]['pakket']),
                    '{}'.format(self.todo[sample]['panel']),
                    '{}'.format(self.todo[sample]['oid']),
                    '{}'.format(PICARD),
                    '{}'.format(GATK),
                    '{}'.format(DBSNP),
                    '{}'.format(REF),
                    '{}'.format(BWA)]

        VARS = ['gDNA', 'cDNA.', 'Eiwitverandering', 'Zygosity',
                'Klassificatie', 'Paraaf staf.', 'Opmerking', 'IGVariant',
                'Paraaf NGS-connaisseur']

        VARS2 = ['Paraaf analist.', 'Gen', 'cDNA.', 'gDNA', 'Eiwitverandering',
                 'Exon', 'Chromosoom', 'NM-nummer (RefSeq)',
                 'Genoom versie', 'Pos.Controle1', 'Pos.Controle2',
                 'Vermelden MLPA kitnr / DL']

        VARS4 = ['Score', 'Over', 'Score', 'Paraaf staf']

        SANGER = ['Gen', 'Chr.', 'g.Start', 'g.Eind', 'Resultaat',
                  'Paraaf analist', 'Paraaf staf']

        MR = MetricsDBReader(self.db_metrics, sample, serie,
                             capture)

        ws1.write(1, 1, sample)
        samenvatting = '''=B4&"("&B2&" + "&B3&"; " &B5&", "&B6&", "&B7&")"'''
        ws1.write(len(MISEQ), 1, samenvatting)

        row1, col1 = write_list(MISEQ, ws1, header='MISEQUENCING',
                                formatheader=headerformat)

        row1, col1 = write_list(SEQRAP, ws1, header='SEQUENCE RAPPORT',
                                formatheader=headerformat, row=row1, col=col1)

        row1, col1 = write_list(MR.get_alignmetrics() + MR.get_hsmetrics(),
                                ws1, row=row1 - len(SEQRAP) - 2, col=1)

        row1, col1 = write_list([200, 30],
                                ws1, row=row1 - 1, col=col1, format=gray)

        row1, col1 = write_list(SNPCHECK, ws1, header='SNPCHECK',
                                formatheader=headerformat, skip=1, row=row1,
                                orientation='cols', format=underlined)

        snpcheck_data = MR.get_snpcheck()
        rs_gpos_dict = get_rs_gpos_dict()
        
        for locus, result in snpcheck_data['COMP'].items():
            rsid = rs_gpos_dict[locus]
            snpcheck_ngs_result = snpcheck_data['NGS'][locus]
            snpcheck_alt_result = snpcheck_data['ALT'][locus]
            if snpcheck_alt_result not in ['WT', 'HET', 'HOM']:
                result = 'NoTaqMan'
                snpcheck_alt_result = 'NoTaqMan'
            snpcheck_excel = [locus, rsid, snpcheck_ngs_result, snpcheck_alt_result, result]
            row1, col1 = write_list(snpcheck_excel,
                                    ws1, row=row1, col=0,
                                    orientation='cols')

        row1, col1 = write_list(INFO, ws1, header='ANALYSE INFO',
                                formatheader=headerformat, row=row1 + 2, col=0)

        row1, col1 = write_list(INFODATA, ws1, row=row1 - len(INFO) - 2, col=1)

        row2, col2 = write_list(SANGER, ws2, formatheader=headerformat,
                                header='NON-CALLABLE SANGER FRAGMENTEN',
                                skip=3, orientation='cols', format=underlined)

        sangers = MR.get_sanger_fragments()
#            df_serie = standaardfragmenten_naar_df(standaardfragmenten)
#            standaard = standaardfragmenten_dict_sample(sample, df_serie)
        standaard = list()

        
        if 'Geen sangers:' in sangers:
            ws2.set_column('A:A', 26)
            sangers = [sangers]
            [sangers.append(' ') for _ in range(len(SANGER))]
            row2, col2 = write_list(sangers, ws2, row=row2-1, col=0,
                                    orientation='cols')

            row2, col2 = write_list([' ', 'Paraaf staf voor gezien: '],
                                    ws2, row=row2, col=0)
        else:
            for line in sangers:
                gene, chrom, start, end, *_ = line
                start = int(start)
                end = int(end)
                score = ' '
                for _ in standaard:
                    if not 'chr.' in _:
                        continue
                    st_chrom = 'chr{}'.format(_['chr.'])
                    st_start = int(_['primer min'])
                    st_end = int(_['primer max'])
                    if st_start <= start and st_end >= end and st_chrom == chrom: 
                        score = _['Score']
                        score = ' '

                out = [gene, chrom, start, end , score]
                row2, col2 = write_list(out, ws2, row=row2, col=0,
                                        orientation='cols')

        ws3.write(3, 11, 'Paraaf', underlined)
        ws3.write(4, 10, 'Eerste beoordeling: ')

        row3, col3 = write_list(VARS, ws3, formatheader=headerformat,
                                header='VARIANTEN IN BRIEF',
                                skip=2, orientation='cols', format=underlined)

        row3, col3 = write_list(VARS2, ws3, formatheader=headerformat,
                                header='GENESIS NOTATIE', row=row3 + 12,
                                orientation='cols', format=underlined)

        for _ in range(12):
            connaisseur_row = int(_) + 4
            connaisseur = '''="chr" & SUBSTITUTE(LOWER(LEFT(A{r},SEARCH("(",A{r},1)-1)), "chr", "") & ":" & MID(A{r}, SEARCH(".",A{r}) + 1,SEARCH(">",A{r}) -  SEARCH(".",A{r}) -2)'''.format(r=connaisseur_row)
            ws3.write(connaisseur_row - 1, 7, connaisseur)
            row = row3 - 13
            gennaam = '''=MID(B{r}, SEARCH("(",B{r},1) + 1,( SEARCH(")",B{r},1) -1 - SEARCH("(",B{r},1) ))'''.format(r=row)
            cnotatie = '''=RIGHT(B{r}, LEN(B{r}) - SEARCH(":",B{r},1))'''.format(r=row)
            gnotatie = '''=RIGHT(A{r}, LEN(A{r}) - SEARCH(":",A{r},1))'''.format(r=row)
            pnotatie = '''=C{r}'''.format(r=row)
            chromosoom = '''=SUBSTITUTE(LOWER(LEFT(A{r},SEARCH("(",A{r},1)-1)), "chr", "")'''.format(r=row)
            nmnummer = '''=LEFT(B{r},SEARCH("(",B{r},1)-1)'''.format(r=row)
            VARS3 = ['', gennaam, cnotatie, gnotatie, pnotatie,  '', chromosoom,
                     nmnummer, 'GRCh37', sample]

            row3, col3 = write_list(VARS3, ws3, formatheader=headerformat,
                                    header=False, row=row3,
                                    orientation='cols')

        row3, col3 = write_list(VARS4, ws3, formatheader=headerformat,
                                header='CNV', row=row3 + 2,
                                orientation='cols', format=underlined)
        ws3.data_validation(row3, 0, row3+3, 0,
                            {'validate': 'list',
                             'source': ['WT', 'NVT', 'AFW', 'NTS']})
        ws3.data_validation(row3, 2, row3+3, 2,
                            {'validate': 'list',
                             'source': ['WT', 'NVT', 'AFW', 'NTS']})

        

        #if self.todo[sample]['cnv_diag']:
        if False:

            ws4 = wb.add_worksheet('CNV')
            ws5 = wb.add_worksheet('CNV info')

            ws4.set_column('A:A', 25)
            ws4.set_column('B:F', 10)
            ws4.set_column('G:G', 11)
            ws4.set_column('H:H', 22)
            ws4.set_column('I:I', 4)
            ws4.set_column('J:J', 19)
            ws4.set_column('K:K', 9)

            ws4.set_column('A:A', 25)
            ws4.set_column('G:G', 18)
            ws4.set_column('I:I', 18)

            cnv_calls = 'output/CNV_{}/Calls/{}.txt'.format(capture, sample)
            cnv_archive = 'output/CNV_{}/archive.txt'.format(capture)
            cnv_excluded = 'output/CNV_{}/excluded.txt'.format(capture)

            if not os.path.isfile(cnv_calls):
                row4, col4 = write_list(['Geen calls', '',
                                         'Paraaf NGS-connaisseur voor gezien: ',
                                         'Paraaf staf voor gezien: '], ws4)
                ws4.set_column('A:A', 35)

            elif os.path.isfile(cnv_calls):
                ws4.write(1, 10, 'Paraaf', underlined)
                ws4.write(2, 9, 'Eerste beoordeling: ')
                with open(cnv_calls, 'r') as f:
                    header = next(f)
                    header = list(header.split())
                    header.insert(0, 'Regio')
                    header.append('Opmerking')
                    header.append('Paraaf NGS-connaisseur')
                    row4, col4 = write_list(header, ws4, formatheader=headerformat,
                                            header='CNV calls', orientation='cols',
                                            format=underlined)

                    for line in f:
                        row4, col4 = write_list(line.split(), ws4, row=row4, orientation='cols')

            with open(cnv_excluded, 'r') as f:
                header = next(f)
                header = list(header.split())
                header.insert(0, 'Regio')
                row5, col5 = write_list(header, ws5, formatheader=headerformat,
                                        header='CNV excluded regions', orientation='cols',
                                        format=underlined)
                for line in f:
                    row5, col5 = write_list(line.split(), ws5, row=row5, orientation='cols')

            with open(cnv_archive, 'r') as f:
                row5, col5 = write_list([], ws5, formatheader=headerformat,
                                        header='CNV archive samples', orientation='cols',
                                        format=underlined, row=row5 + 4)
                for line in f:
                    row5, col5 = write_list(line.split(), ws5, row=row5, orientation='cols')
            
        
        
        
        

    