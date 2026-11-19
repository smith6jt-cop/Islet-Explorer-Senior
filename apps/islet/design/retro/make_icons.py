import math, sys
sys.path.insert(0, sys.argv[1] + "/logo")
from make_logo import *
OUT = f"{S}/logo/retro"
def svg(body, w=512, h=512, label="", fonts=False):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'role="img" aria-label="{label}">{style() if fonts else ""}{body}</svg>')
def plate(fill=CREAM):
    return (f'<circle cx="256" cy="256" r="250" fill="{INK}"/><circle cx="256" cy="256" r="238" fill="{MUST}"/>'
            f'<circle cx="256" cy="256" r="222" fill="{fill}"/>')
def blob(R, n=72, k=7, amp=.035, ph=0.0):
    pts = []
    for i in range(n):
        a = 2*math.pi*i/n; r = R*(1+amp*math.sin(k*a+ph))
        pts.append(f"{r*math.cos(a):.1f},{r*math.sin(a):.1f}")
    return " ".join(pts)
def sunrays(cx, cy, r, n=20, fill="#efe0b4"):
    d = []
    for i in range(n):
        a0 = 2*math.pi*i/n; a1 = 2*math.pi*(i+.5)/n
        d.append(f"M{cx},{cy} L{cx+r*math.cos(a0):.1f},{cy+r*math.sin(a0):.1f} L{cx+r*math.cos(a1):.1f},{cy+r*math.sin(a1):.1f}Z")
    return f'<g clip-path="url(#pl)"><path fill="{fill}" d="{" ".join(d)}"/></g>'
CLIP = '<defs><clipPath id="pl"><circle cx="256" cy="256" r="222"/></clipPath></defs>'
I = {}
# --- states (no character) ---------------------------------------------------------------
I["state-welcome"] = svg(CLIP + plate() + sunrays(256, 215, 300)
    + f'<g transform="translate(-74 -5) scale(.55)">{hero()}</g>'
    + f'<path fill="{REDD}" d="M64,356 L120,356 L120,410 L64,410 L84,383Z M448,356 L392,356 L392,410 L448,410 L428,383Z"/>'
    + f'<rect x="100" y="338" width="312" height="62" fill="{RED}"/>'
    + f'<text class="c" x="256" y="388" font-size="52" letter-spacing="5" text-anchor="middle" fill="{CREAM}">WELCOME!</text>',
    label="Welcome", fonts=True)
dots = "".join(f'<circle cx="{x}" cy="{y}" r="5" fill="{INK}" fill-opacity=".3"/>' for x in range(110, 420, 36) for y in range(120, 420, 36))
I["state-nothing-selected"] = svg(CLIP + plate() + f'<g clip-path="url(#pl)">{dots}</g>'
    + f'<path d="M296,300 L384,396" stroke="{INK}" stroke-width="30" stroke-linecap="round"/>'
    + f'<path d="M318,324 L384,396" stroke="{RED}" stroke-width="16" stroke-linecap="round"/>'
    + f'<circle cx="226" cy="228" r="104" fill="{WHITE}" fill-opacity=".85" stroke="{INK}" stroke-width="16"/>'
    + f'<circle cx="226" cy="228" r="52" fill="none" stroke="{INK}" stroke-opacity=".55" stroke-width="7" stroke-dasharray="14 12"/>'
    + f'<path d="M158,196 Q172,156 212,146" fill="none" stroke="{TURQ}" stroke-width="10" stroke-linecap="round"/>'
    + star(404, 150, 18, MUST), label="Nothing selected yet")
I["state-unavailable"] = svg(CLIP + plate()
    + f'<circle cx="256" cy="262" r="140" fill="{WHITE}" fill-opacity=".6" stroke="{INK}" stroke-width="9" stroke-dasharray="22 16" stroke-linecap="round"/>'
    + f'<text class="c" x="256" y="340" font-size="220" text-anchor="middle" fill="{RED}">?</text>'
    + star(110, 120, 16, MUST) + star(402, 120, 16, MUST), label="Unavailable — not a zero", fonts=True)
I["state-working"] = svg(CLIP + plate() + sunrays(256, 256, 300)
    + f'<path d="M186,132 H326 Q326,212 270,256 Q326,300 326,388 H186 Q186,300 242,256 Q186,212 186,132Z" fill="{WHITE}" stroke="{INK}" stroke-width="9" stroke-linejoin="round"/>'
    + f'<path d="M208,176 H304 Q296,216 262,244 H250 Q216,216 208,176Z" fill="{MAG}"/>'
    + f'<path d="M256,244 V372" stroke="{MAG}" stroke-width="6" stroke-dasharray="10 8"/>'
    + f'<path d="M204,384 Q256,296 308,384Z" fill="{ORA}"/>'
    + f'<path d="M160,122 H352 M160,398 H352" stroke="{INK}" stroke-width="24" stroke-linecap="round"/>'
    + f'<path d="M160,122 H352 M160,398 H352" stroke="{MUST}" stroke-width="10" stroke-linecap="round"/>'
    + f'<path d="M376,200 Q400,256 376,312 M136,200 Q112,256 136,312" fill="none" stroke="{RED}" stroke-width="7" stroke-linecap="round"/>',
    label="Working on it")
I["state-whoops"] = svg(CLIP + plate()
    + f'<path d="M150,388 Q170,350 240,356 Q300,330 350,356 Q420,360 408,392 Q380,420 280,414 Q170,424 150,388Z" fill="{MAG}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>'
    + f'<circle cx="126" cy="372" r="10" fill="{MAG}" stroke="{INK}" stroke-width="6"/><circle cx="430" cy="352" r="8" fill="{MAG}" stroke="{INK}" stroke-width="6"/>'
    + f'<g transform="translate(206 300) rotate(-66)"><path d="M-18,-112 h36 v54 l52,112 q8,20 -14,20 h-112 q-22,0 -14,-20 l52,-112Z" fill="{WHITE}" stroke="{INK}" stroke-width="9" stroke-linejoin="round"/>'
    + f'<path d="M-52,20 h104 l18,36 q6,14 -10,14 h-120 q-16,0 -10,-14Z" fill="{MAG}"/><path d="M-26,-112 h52" stroke="{INK}" stroke-width="12" stroke-linecap="round"/></g>'
    + f'<g transform="rotate(10 366 170)">{seal(366, 170, 84, 12, MUST)}<text class="c" x="366" y="212" font-size="120" text-anchor="middle" fill="{RED}">!</text></g>',
    label="Whoops — something went wrong", fonts=True)
# --- please stand by ---------------------------------------------------------------------
tv = f'''<path d="M196,104 L256,150 L330,92" fill="none" stroke="{INK}" stroke-width="8" stroke-linecap="round"/>
<circle cx="196" cy="104" r="9" fill="{RED}"/><circle cx="330" cy="92" r="9" fill="{RED}"/>
<path d="M150,392 L128,436 M362,392 L384,436" stroke="{INK}" stroke-width="12" stroke-linecap="round"/>
<rect x="88" y="148" width="336" height="250" rx="34" fill="{RED}" stroke="{INK}" stroke-width="8"/>
<rect x="112" y="172" width="230" height="202" rx="40" fill="{TEAL}" stroke="{INK}" stroke-width="8"/>
<circle cx="227" cy="273" r="78" fill="none" stroke="{CREAM}" stroke-width="6"/>
<circle cx="227" cy="273" r="50" fill="none" stroke="{MUST}" stroke-width="6"/>
<path d="M132,273 L322,273 M227,186 L227,360" stroke="{CREAM}" stroke-width="4"/>
<rect x="126" y="252" width="202" height="44" fill="{INK}"/>
<text class="c" x="227" y="286" font-size="34" letter-spacing="2" text-anchor="middle" fill="{CREAM}">PLEASE STAND BY</text>
<circle cx="382" cy="212" r="20" fill="{CREAM}" stroke="{INK}" stroke-width="6"/><path d="M382,212 L382,196" stroke="{INK}" stroke-width="6" stroke-linecap="round"/>
<circle cx="382" cy="270" r="20" fill="{CREAM}" stroke="{INK}" stroke-width="6"/><path d="M382,270 L394,260" stroke="{INK}" stroke-width="6" stroke-linecap="round"/>
<path d="M366,322 H398 M366,338 H398 M366,354 H398" stroke="{INK}" stroke-width="5" stroke-linecap="round"/>'''
I["please-stand-by"] = svg(CLIP + plate() + sunrays(256, 273, 300) + tv, label="Please stand by", fonts=True)
# --- tab icons ---------------------------------------------------------------------------
def wedge(cx, cy, r, a0, a1, fill, dx=0, dy=0):
    x0, y0 = cx+r*math.cos(a0), cy+r*math.sin(a0); x1, y1 = cx+r*math.cos(a1), cy+r*math.sin(a1)
    big = 1 if a1-a0 > math.pi else 0
    return (f'<path transform="translate({dx} {dy})" d="M{cx},{cy} L{x0:.1f},{y0:.1f} A{r},{r} 0 {big} 1 {x1:.1f},{y1:.1f}Z" '
            f'fill="{fill}" stroke="{INK}" stroke-width="7" stroke-linejoin="round"/>')
d2r = math.radians
pie = (f'<ellipse cx="250" cy="290" rx="176" ry="176" fill="{WHITE}" stroke="{INK}" stroke-width="7"/>'
       f'<polygon transform="translate(250 290)" points="{blob(150, 96, 22, .03)}" fill="{TAN}" stroke="{INK}" stroke-width="7"/>'
       + wedge(250, 290, 128, d2r(-40), d2r(80), ORA) + wedge(250, 290, 128, d2r(80), d2r(170), TURQ)
       + wedge(250, 290, 128, d2r(170), d2r(250), MUST) + wedge(250, 290, 128, d2r(250), d2r(320), MAG, 22, -30)
       + f'<path d="M190,120 q-18,-22 0,-44 q18,-22 0,-44 M250,104 q-18,-22 0,-44 q18,-22 0,-40" fill="none" stroke="{INK}" stroke-opacity=".5" stroke-width="7" stroke-linecap="round"/>')
I["tab-composition"] = svg(CLIP + plate() + pie, label="Composition — who gets what slice")
rocket = f'''<path d="M92,404 Q150,250 300,214" fill="none" stroke="{INK}" stroke-width="7" stroke-dasharray="2 20" stroke-linecap="round"/>
<circle cx="92" cy="404" r="14" fill="{MAG}" stroke="{INK}" stroke-width="6"/><circle cx="142" cy="318" r="12" fill="{ORA}" stroke="{INK}" stroke-width="6"/>
<circle cx="214" cy="252" r="10" fill="{MUST}" stroke="{INK}" stroke-width="6"/>
<g transform="translate(340 180) rotate(52)">
<path d="M-22,74 Q0,140 22,74Z" fill="{ORA}" stroke="{INK}" stroke-width="6" stroke-linejoin="round"/>
<path d="M-40,30 L-78,92 L-34,74Z M40,30 L78,92 L34,74Z" fill="{RED}" stroke="{INK}" stroke-width="7" stroke-linejoin="round"/>
<path d="M0,-110 Q56,-40 38,78 L-38,78 Q-56,-40 0,-110Z" fill="{WHITE}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>
<path d="M0,-110 Q26,-78 36,-50 L-36,-50 Q-26,-78 0,-110Z" fill="{RED}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>
<circle cx="0" cy="-4" r="22" fill="{TURQ}" stroke="{INK}" stroke-width="8"/></g>'''
I["tab-trajectory"] = svg(CLIP + plate() + rocket + star(120, 150, 18, MUST) + star(400, 340, 22, MUST) + star(330, 420, 12, RED),
                          label="Trajectory — where the cells are headed")
g = lambda x: math.exp(-((x-256)/62)**2/2)
curve = " ".join(f"L{x},{388-190*g(x):.1f}" for x in range(86, 427, 10))
bell = f'''<path d="M86,388 {curve} L426,388Z" fill="{TURQ}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>
<path d="M70,388 H442" stroke="{INK}" stroke-width="8" stroke-linecap="round"/>
<path d="M256,388 V214" stroke="{INK}" stroke-width="5" stroke-dasharray="8 10"/>
<g transform="translate(256 150)"><path d="M0,-52 V-70" stroke="{INK}" stroke-width="8" stroke-linecap="round"/>
<path d="M-50,34 Q-40,20 -38,-10 Q-36,-52 0,-52 Q36,-52 38,-10 Q40,20 50,34Z" fill="{MUST}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>
<circle cx="0" cy="44" r="11" fill="{RED}" stroke="{INK}" stroke-width="6"/>
<path d="M-74,-30 Q-86,-6 -74,18 M74,-30 Q86,-6 74,18 M-96,-40 Q-114,-6 -96,28 M96,-40 Q114,-6 96,28" fill="none" stroke="{RED}" stroke-width="6" stroke-linecap="round"/></g>'''
I["tab-statistics"] = svg(CLIP + plate() + bell, label="Statistics — the bell curve, with bell")
mp = f'''<path d="M96,170 L196,140 L316,176 L416,146 L416,372 L316,402 L196,366 L96,396Z" fill="{WHITE}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>
<path d="M196,140 V366 M316,176 V402" stroke="{INK}" stroke-opacity=".35" stroke-width="5"/>
<path d="M130,350 Q170,270 230,300 T330,250" fill="none" stroke="{RED}" stroke-width="7" stroke-dasharray="3 16" stroke-linecap="round"/>
<path d="M124,338 l14,24 M138,338 l-14,24" stroke="{RED}" stroke-width="7" stroke-linecap="round"/>
<path d="M340,262 Q286,190 286,150 A54,54 0 1 1 394,150 Q394,190 340,262Z" fill="{RED}" stroke="{INK}" stroke-width="8" stroke-linejoin="round"/>
<circle cx="340" cy="150" r="22" fill="{CREAM}" stroke="{INK}" stroke-width="7"/>'''
I["tab-spatial"] = svg(CLIP + plate() + mp + star(150, 210, 14, MUST), label="Spatial — you are here")
def house(x, y, w, c, roof):
    return (f'<rect x="{x}" y="{y}" width="{w}" height="{w*.7:.0f}" fill="{c}" stroke="{INK}" stroke-width="7"/>'
            f'<path d="M{x-14},{y} L{x+w/2},{y-w*.5:.0f} L{x+w+14},{y}Z" fill="{roof}" stroke="{INK}" stroke-width="7" stroke-linejoin="round"/>'
            f'<rect x="{x+w*.38:.0f}" y="{y+w*.3:.0f}" width="{w*.24:.0f}" height="{w*.4:.0f}" fill="{INK}"/>'
            f'<rect x="{x+w*.1:.0f}" y="{y+w*.18:.0f}" width="{w*.2:.0f}" height="{w*.2:.0f}" fill="{WHITE}" stroke="{INK}" stroke-width="5"/>')
fence = "".join(f'<path d="M{x},372 V338 l8,-10 l8,10 V372Z" fill="{WHITE}" stroke="{INK}" stroke-width="4" stroke-linejoin="round"/>' for x in range(70, 440, 22))
hood = (f'<circle cx="372" cy="132" r="34" fill="{MUST}" stroke="{INK}" stroke-width="7"/>'
        + house(84, 250, 104, TURQ, RED) + house(206, 232, 116, MUST, TEAL) + house(340, 252, 100, MAG, ORA)
        + f'<rect x="40" y="366" width="432" height="60" fill="#8cc9a0"/>' + fence)
I["tab-neighborhoods"] = svg(CLIP + plate() + f'<g clip-path="url(#pl)">{hood}</g>', label="Neighborhoods — meet the neighbours")
venn = f'''<circle cx="196" cy="262" r="118" fill="{MAG}" fill-opacity=".9" stroke="{INK}" stroke-width="8"/>
<circle cx="316" cy="262" r="118" fill="{TURQ}" fill-opacity=".85" stroke="{INK}" stroke-width="8"/>
<path d="M256,300 q-26,-22 -14,-38 q10,-10 14,4 q4,-14 14,-4 q12,16 -14,38Z" fill="{RED}" stroke="{INK}" stroke-width="5" stroke-linejoin="round"/>'''
I["tab-crossmodal"] = svg(CLIP + plate() + venn + star(256, 120, 22, MUST), label="Cross-modal — a perfect match")
scope = f'''<path d="M150,410 H370" stroke="{INK}" stroke-width="26" stroke-linecap="round"/>
<path d="M300,400 Q392,330 330,214" fill="none" stroke="{INK}" stroke-width="24" stroke-linecap="round"/>
<path d="M176,356 H300" stroke="{INK}" stroke-width="14" stroke-linecap="round"/>
<g transform="rotate(-28 262 220)"><rect x="232" y="104" width="60" height="170" rx="12" fill="{RED}" stroke="{INK}" stroke-width="8"/>
<rect x="240" y="78" width="44" height="34" rx="6" fill="{CREAM}" stroke="{INK}" stroke-width="8"/>
<rect x="244" y="274" width="36" height="34" fill="{CREAM}" stroke="{INK}" stroke-width="8"/></g>
<circle cx="340" cy="250" r="20" fill="{MUST}" stroke="{INK}" stroke-width="7"/>'''
I["tab-viewer"] = svg(CLIP + plate() + scope + star(206, 318, 26, MUST) + star(140, 190, 14, RED), label="Viewer — take a closer look")
book = f'''<g transform="rotate(-8 256 256)"><rect x="150" y="96" width="232" height="320" rx="8" fill="{INK}"/>
<rect x="138" y="88" width="232" height="320" rx="8" fill="{TEAL}" stroke="{INK}" stroke-width="8"/>
<rect x="158" y="108" width="192" height="280" fill="none" stroke="{CREAM}" stroke-width="4"/>
<text class="c" x="254" y="290" font-size="190" text-anchor="middle" fill="{MUST}">?</text>
<text class="c" x="254" y="356" font-size="50" letter-spacing="4" text-anchor="middle" fill="{CREAM}">HELP</text></g>'''
I["help-booklet"] = svg(CLIP + plate() + book + star(410, 150, 20, RED), label="Help", fonts=True)
# --- stamps and ornaments ----------------------------------------------------------------
def stamp(text, colour):
    w = 110 + 44*len(text)
    return svg(f'''<g transform="rotate(-7 {w/2} 120)" fill="none" stroke="{colour}">
<rect x="26" y="48" width="{w-52}" height="144" rx="18" stroke-width="12"/>
<rect x="44" y="66" width="{w-88}" height="108" rx="8" stroke-width="4"/>
<text class="c" x="{w/2}" y="156" font-size="104" letter-spacing="6" text-anchor="middle" fill="{colour}" stroke="none">{text}</text></g>''',
               w, 240, text.title(), fonts=True)
I["stamp-provisional"] = stamp("PROVISIONAL", RED)
I["stamp-withheld"] = stamp("WITHHELD", RED)
I["stamp-not-a-zero"] = stamp("NOT A ZERO", "#c98a12")
I["stamp-approved"] = stamp("APPROVED", "#1a8f88")
def burst8(cx, cy, r, c):
    return "".join(f'<path d="M{cx},{cy-r} L{cx},{cy+r}" stroke="{c}" stroke-width="5" stroke-linecap="round" transform="rotate({a} {cx} {cy})"/>' for a in (0, 45, 90, 135)) \
        + f'<circle cx="{cx}" cy="{cy}" r="7" fill="{c}"/>'
I["divider-atomic"] = svg(
    f'<path d="M20,30 H520 M680,30 H1180" stroke="{MUST}" stroke-width="4" stroke-linecap="round"/>'
    f'<path d="M60,42 H500 M700,42 H1140" stroke="{TURQ}" stroke-width="2" stroke-linecap="round"/>'
    + burst8(600, 30, 24, RED) + star(545, 30, 14, MUST) + star(655, 30, 14, MUST), 1200, 60, "")
for k, v in I.items():
    open(f"{OUT}/{k}.svg", "w").write(v)
print(len(I))
