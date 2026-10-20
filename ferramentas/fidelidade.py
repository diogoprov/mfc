import re, sys, unicodedata
def norm(s):
    s = unicodedata.normalize("NFKD", s.lower()); return "".join(c for c in s if not unicodedata.combining(c))
def slides(path):
    out = {}; cur = None
    for line in open(path):
        m = re.match(r"=== Slide (\d+)", line)
        if m: cur = int(m.group(1)); out[cur] = []; continue
        if cur and not line.strip().startswith("[") : out[cur].append(line.strip())
    return out
stop = set("porque quando sobre entre outras outros também apenas através enquanto durante depois antes muito mesmo ainda podem pode deve devem sendo esses essas estes estas aquele aquela dados atributo atributos espécies modelo modelos maneira caracter caracteres trait traits figure slide".split())
orig, ini, fim, deck = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
D = norm(open(deck).read())
for n, linhas in slides(orig).items():
    if not (ini <= n <= fim): continue
    txt = norm(" ".join(linhas))
    pal = [w for w in set(re.findall(r"[a-z]{7,}", txt)) if w not in stop]
    if len(pal) < 4: continue
    falta = [w for w in pal if w[:6] not in D]
    cob = 1 - len(falta) / len(pal)
    if cob < .6: print(f"slide {n}: {cob:.0%} | faltam: {', '.join(sorted(falta)[:12])}")
