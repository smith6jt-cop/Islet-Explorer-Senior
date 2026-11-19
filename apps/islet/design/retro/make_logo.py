import base64, io, math, sys
from fontTools import subset
from fontTools.ttLib import TTFont
S = sys.argv[1]
def b64font(path, text):
    f = TTFont(path)
    opt = subset.Options(); opt.layout_features = ['*']
    s = subset.Subsetter(opt); s.populate(text=text); s.subset(f)
    buf = io.BytesIO(); f.flavor = 'woff2' if False else None; f.save(buf)
    return base64.b64encode(buf.getvalue()).decode()
SCRIPT = b64font(f"{S}/fonts/Lobster-Regular.ttf", "MultiPlex")
CAPS = b64font(f"{S}/fonts/BebasNeue-Regular.ttf", "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 -!,.&'")
TEAL, TEAL2, TEALD = "#0d5461", "#136676", "#083a44"
CREAM, MUST, RED, REDD = "#f6ecd0", "#f0b93a", "#d8452b", "#96301d"
MAG, ORA, TURQ = "#d8246e", "#f58a23", "#1fb9b0"
def star(x, y, r, fill=CREAM):
    return (f'<path fill="{fill}" d="M{x},{y-r} Q{x},{y} {x+r},{y} Q{x},{y} {x},{y+r} '
            f'Q{x},{y} {x-r},{y} Q{x},{y} {x},{y-r}Z"/>')
def style():
    return f'''<style>
@font-face{{font-family:"LogoScript";src:url(data:font/ttf;base64,{SCRIPT}) format("truetype");}}
@font-face{{font-family:"LogoCaps";src:url(data:font/ttf;base64,{CAPS}) format("truetype");}}
.s{{font-family:"LogoScript";}} .c{{font-family:"LogoCaps";}}
</style>'''
def burst(cx, cy, r, n=24):
    out = []
    for i in range(n):
        a0 = 2*math.pi*(i/n); a1 = 2*math.pi*((i+.5)/n)
        out.append(f'M{cx},{cy} L{cx+r*math.cos(a0):.1f},{cy+r*math.sin(a0):.1f} '
                   f'L{cx+r*math.cos(a1):.1f},{cy+r*math.sin(a1):.1f}Z')
    return f'<path fill="{TEAL2}" d="{" ".join(out)}"/>'
INK, PINK, TAN, WHITE = TEALD, "#f29a8e", "#d9a85c", "#fffaf0"

def cell(x, y, r, fill):
    return (f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{CREAM}" stroke-width="6"/>'
            f'<circle cx="{x+r*.18:.1f}" cy="{y+r*.12:.1f}" r="{r*.42:.1f}" fill="{INK}" fill-opacity=".55"/>'
            f'<path d="M{x-r*.55:.1f},{y-r*.3:.1f} Q{x-r*.4:.1f},{y-r*.62:.1f} {x-r*.05:.1f},{y-r*.66:.1f}" fill="none" '
            f'stroke="{WHITE}" stroke-opacity=".7" stroke-width="{r*.14:.1f}" stroke-linecap="round"/>')
def hero():
    """The main figure: a cluster of cells in every colour, with more drifting around it."""
    o = []
    for x, y, r, c in ((318, 318, 40, MAG), (296, 452, 30, MUST), (392, 540, 24, ORA), (372, 222, 26, TURQ),
                       (884, 300, 38, ORA), (902, 436, 30, MAG), (812, 544, 26, MUST), (826, 216, 24, TURQ)):
        o.append(cell(x, y, r, c))
    ring = [MAG, ORA, MUST, TURQ, MAG, ORA, MUST]
    for i, c in enumerate(ring):
        a = 2*math.pi*i/len(ring) - math.pi/2
        o.append(cell(600+124*math.cos(a), 400+124*math.sin(a), 52, c))
    o.append(cell(600, 400, 66, TURQ))
    return "\n".join(o)
def wordmark(dy=0):
    return f'''<g transform="translate(0 {dy})">
<g transform="rotate(-5 600 640)">
<text class="s" x="606" y="668" font-size="200" text-anchor="middle" fill="{REDD}">MultiPlex</text>
<text class="s" x="600" y="660" font-size="200" text-anchor="middle" fill="{CREAM}">MultiPlex</text>
</g>
<g transform="translate(600 862) skewX(-9)">
<text class="c" x="9" y="9" font-size="236" letter-spacing="14" text-anchor="middle" fill="{TEALD}">5000</text>
<text class="c" x="0" y="0" font-size="236" letter-spacing="14" text-anchor="middle" fill="{MUST}">5000</text>
</g></g>'''
def full():
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 1200" width="1200" height="1200" role="img" aria-label="MultiPlex5000">
{style()}
<defs>
<clipPath id="field"><circle cx="600" cy="600" r="462"/></clipPath>
<path id="arcTop" d="M103,600 A497,497 0 0 1 1097,600"/>
<path id="arcBot" d="M64,600 A536,536 0 0 0 1136,600"/>
</defs>
<circle cx="600" cy="600" r="592" fill="{TEALD}"/>
<circle cx="600" cy="600" r="578" fill="{CREAM}"/>
<circle cx="600" cy="600" r="564" fill="{MUST}"/>
<circle cx="600" cy="600" r="472" fill="{CREAM}"/>
<circle cx="600" cy="600" r="462" fill="{TEAL}"/>
<g clip-path="url(#field)">{burst(600, 400, 900)}</g>
<text class="c" font-size="70" letter-spacing="9" fill="{TEALD}" text-anchor="middle"><textPath href="#arcTop" startOffset="50%">PROTEINS, GENES &amp; MORE</textPath></text>
{star(600, 1116, 24, TEALD)}{star(478, 1101, 16, TEALD)}{star(722, 1101, 16, TEALD)}
<g transform="translate(120 706) rotate(-14) scale(.8) translate(-742 -664)">{seal(742, 664, 92)}<circle cx="742" cy="664" r="62" fill="none" stroke="{CREAM}" stroke-width="3"/><text class="c" x="742" y="656" font-size="46" letter-spacing="2" text-anchor="middle" fill="{CREAM}">ALL</text><text class="c" x="742" y="700" font-size="46" letter-spacing="2" text-anchor="middle" fill="{CREAM}">NEW!</text></g>{star(1072, 720, 30, TEALD)}
<g transform="translate(600 366) scale(.88) translate(-600 -400)">{hero()}</g>
{star(300, 250, 26)}{star(250, 600, 18, MUST)}{star(985, 600, 20, MUST)}{star(860, 215, 16)}
{wordmark(8)}
<g>
<path fill="{REDD}" d="M96,930 L212,930 L212,1012 L96,1012 L128,971Z"/>
<path fill="{REDD}" d="M1104,930 L988,930 L988,1012 L1104,1012 L1072,971Z"/>
<path fill="#5f1d11" d="M172,980 L212,1012 L212,980Z M1028,980 L988,1012 L988,980Z"/>
<rect x="172" y="896" width="856" height="84" fill="{RED}"/>
<text class="c" x="600" y="961" font-size="64" letter-spacing="5" text-anchor="middle" fill="{CREAM}">TOMORROW'S SCIENCE, TODAY!</text>
</g>
</svg>'''
def mark():
    # compact form for the navigation rail: no small lettering
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 1200" width="1200" height="1200" role="img" aria-label="MultiPlex5000">
{style()}
<defs><clipPath id="field"><circle cx="600" cy="600" r="540"/></clipPath></defs>
<circle cx="600" cy="600" r="592" fill="{TEALD}"/>
<circle cx="600" cy="600" r="578" fill="{MUST}"/>
<circle cx="600" cy="600" r="552" fill="{CREAM}"/>
<circle cx="600" cy="600" r="540" fill="{TEAL}"/>
<g clip-path="url(#field)">{burst(600, 430, 900)}</g>
<g transform="translate(-66 -30) scale(1.11)">{hero()}</g>
<g transform="translate(-90 20) scale(1.15)">{wordmark(0)}</g>
</svg>'''

def seal(cx, cy, r, n=16, fill=RED):
    pts = []
    for i in range(2*n):
        rr = r if i % 2 == 0 else r*0.8
        a = math.pi*i/n
        pts.append(f"{cx+rr*math.cos(a):.1f},{cy+rr*math.sin(a):.1f}")
    return f'<polygon fill="{fill}" points="{" ".join(pts)}"/>'
def pamphlet():
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 1200" width="900" height="1200" role="img" aria-label="MultiPlex5000">
{style()}
<defs><clipPath id="panel"><rect x="56" y="56" width="788" height="600"/></clipPath></defs>
<rect width="900" height="1200" rx="26" fill="{TEALD}"/>
<rect x="18" y="18" width="864" height="1164" rx="14" fill="{CREAM}"/>
<rect x="38" y="38" width="824" height="1124" fill="none" stroke="{TEALD}" stroke-width="5"/>
<rect x="56" y="56" width="788" height="600" fill="{TEAL}"/>
<g clip-path="url(#panel)">{burst(450, 400, 1000)}
<g transform="translate(-140 10)">{hero()}</g>
{star(130, 250, 22)}{star(770, 180, 18, MUST)}{star(700, 610, 14)}
</g>
<text class="c" x="450" y="128" font-size="46" letter-spacing="8" text-anchor="middle" fill="{CREAM}">A HANDY GUIDE TO YOUR TISSUE</text>
<rect x="56" y="656" width="788" height="14" fill="{MUST}"/>
<g transform="rotate(-4 450 800)">
<text class="s" x="446" y="836" font-size="196" text-anchor="middle" fill="{RED}">MultiPlex</text>
<text class="s" x="440" y="829" font-size="196" text-anchor="middle" fill="{TEAL}">MultiPlex</text>
</g>
<g transform="translate(462 1022) skewX(-9)">
<text class="c" x="7" y="7" font-size="212" letter-spacing="12" text-anchor="middle" fill="{TEALD}">5000</text>
<text class="c" font-size="212" letter-spacing="12" text-anchor="middle" fill="{RED}">5000</text>
</g>
<rect x="56" y="1062" width="788" height="82" fill="{MUST}"/>
<text class="c" x="450" y="1122" font-size="50" letter-spacing="5" text-anchor="middle" fill="{TEALD}">PROTEINS, GENES &amp; MORE</text>
<g transform="translate(152 572) rotate(-12) scale(.86) translate(-742 -664)">{seal(742, 664, 92)}
<circle cx="742" cy="664" r="62" fill="none" stroke="{CREAM}" stroke-width="3"/>
<text class="c" x="742" y="656" font-size="46" letter-spacing="2" text-anchor="middle" fill="{CREAM}">ALL</text>
<text class="c" x="742" y="700" font-size="46" letter-spacing="2" text-anchor="middle" fill="{CREAM}">NEW!</text></g>
</svg>"""
def nameplate():
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1500 800" width="1500" height="800" role="img" aria-label="MultiPlex5000">
{style()}
<defs><clipPath id="disc"><circle cx="340" cy="400" r="270"/></clipPath></defs>
<rect width="1500" height="800" rx="40" fill="{TEALD}"/>
<rect x="16" y="16" width="1468" height="768" rx="28" fill="{CREAM}"/>
<rect x="36" y="36" width="1428" height="728" rx="14" fill="none" stroke="{TEALD}" stroke-width="4"/>
<circle cx="340" cy="400" r="296" fill="{TEALD}"/>
<circle cx="340" cy="400" r="284" fill="{MUST}"/>
<circle cx="340" cy="400" r="270" fill="{TEAL}"/>
<g clip-path="url(#disc)">{burst(340, 400, 400)}</g>
<g transform="translate(-72 120) scale(.70)">{hero()}</g>
<g transform="rotate(-4 1040 330)">
<text class="s" x="666" y="388" font-size="212" fill="{RED}">MultiPlex</text>
<text class="s" x="660" y="380" font-size="212" fill="{TEAL}">MultiPlex</text>
</g>
<g transform="translate(672 632) skewX(-9)">
<text class="c" x="8" y="8" font-size="236" letter-spacing="12" fill="{TEALD}">5000</text>
<text class="c" font-size="236" letter-spacing="12" fill="{RED}">5000</text>
</g>
<rect x="1128" y="462" width="8" height="172" fill="{MUST}"/>
<text class="c" font-size="66" letter-spacing="6" fill="{TEAL}">
<tspan x="1156" y="510">PROTEINS,</tspan><tspan x="1156" y="571">GENES</tspan><tspan x="1156" y="632">&amp; MORE</tspan></text>
{star(690, 716, 20, MUST)}{star(1410, 716, 20, MUST)}
<text class="c" x="1050" y="734" font-size="52" letter-spacing="7" text-anchor="middle" fill="{TEALD}">TOMORROW'S SCIENCE, TODAY!</text>
<g transform="translate(1368 136) rotate(12) scale(.92) translate(-742 -664)">{seal(742, 664, 92)}<circle cx="742" cy="664" r="62" fill="none" stroke="{CREAM}" stroke-width="3"/><text class="c" x="742" y="656" font-size="46" letter-spacing="2" text-anchor="middle" fill="{CREAM}">ALL</text><text class="c" x="742" y="700" font-size="46" letter-spacing="2" text-anchor="middle" fill="{CREAM}">NEW!</text></g>
</svg>"""
for name, fn in (("logo", full), ("mark", mark), ("pamphlet", pamphlet), ("nameplate", nameplate)):
    open(f"{S}/logo/multiplex5000-{name}.svg", "w").write(fn())
