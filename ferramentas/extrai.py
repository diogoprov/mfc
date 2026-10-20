import sys, re, os, zipfile, xml.etree.ElementTree as ET
src, out = sys.argv[1], sys.argv[2]
z = zipfile.ZipFile(src)
names = z.namelist()
A = '{http://schemas.openxmlformats.org/drawingml/2006/main}'
def texts(xml):
    root = ET.fromstring(xml)
    paras = []
    for p in root.iter(A+'p'):
        t = ''.join(x.text or '' for x in p.iter(A+'t')).strip()
        if t: paras.append(t)
    return paras
# slide order from presentation.xml rels
pres = z.read('ppt/presentation.xml').decode()
rels = z.read('ppt/_rels/presentation.xml.rels').decode()
rmap = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]+)"', rels))
rmap.update({b:a for a,b in re.findall(r'Target="([^"]+)"[^>]*Id="(rId\d+)"', rels)})
ids = re.findall(r'<p:sldId [^>]*r:id="(rId\d+)"', pres)
order = [rmap[i] for i in ids]
with open(out,'w') as f:
    for k, tgt in enumerate(order, 1):
        sp = 'ppt/' + tgt
        f.write(f'\n=== Slide {k} ({os.path.basename(sp)})\n')
        for t in texts(z.read(sp)): f.write('  ' + t + '\n')
        rp = sp.replace('slides/', 'slides/_rels/') + '.rels'
        if rp in names:
            r = z.read(rp).decode()
            media = re.findall(r'Target="\.\./media/([^"]+)"', r)
            notes = re.findall(r'Target="\.\./notesSlides/([^"]+)"', r)
            links = re.findall(r'Target="(https?://[^"]+)"', r)
            if media: f.write('  [mídia] ' + ', '.join(media) + '\n')
            if links: f.write('  [links] ' + ' '.join(links) + '\n')
            for n in notes:
                nt = [t for t in texts(z.read('ppt/notesSlides/'+n)) if not re.fullmatch(r'\d+', t)]
                if nt: f.write('  [notas] ' + ' / '.join(nt) + '\n')
print(out, len(order))
