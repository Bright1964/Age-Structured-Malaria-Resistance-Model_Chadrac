"""Static checks only. Does not execute R or validate epidemiology."""
from pathlib import Path
import re,csv
root=Path(__file__).resolve().parents[1]
scripts=sorted((root/'scripts').glob('*.R'))

def brackets(text):
    pairs={')':'(', ']':'[', '}':'{'}
    stack=[];string=None;escape=False;comment=False;ln=1;col=0
    for ch in text:
        if ch=='\n': ln+=1;col=0;comment=False
        else: col+=1
        if comment:continue
        if string:
            if escape:escape=False
            elif ch=='\\':escape=True
            elif ch==string:string=None
            continue
        if ch=='#':comment=True;continue
        if ch in '"\'`':string=ch;continue
        if ch in '([{':stack.append((ch,ln,col))
        elif ch in ')]}':
            if not stack or stack[-1][0]!=pairs[ch]:
                raise ValueError(f'mismatched {ch} at {ln}:{col}; previous={stack[-1] if stack else None}')
            stack.pop()
    if string:raise ValueError(f'unclosed {string} string')
    if stack:raise ValueError(f'unclosed bracket {stack[-1]}')

for p in [root/'run_all.R',root/'install_packages.R',*scripts]:
    brackets(p.read_text(encoding='utf-8'))
    print('balanced:',p.relative_to(root))

main=(root/'scripts/01_simulate_main.R').read_text()
province=(root/'scripts/04_figure2_provinces.R').read_text()
figs=(root/'scripts/02_figures_main.R').read_text()
supp=(root/'scripts/03_figureS1.R').read_text()
assert 'DT <- 1.0' in main and 'DT <- 1.0' in province
assert 'THRESHOLD_COVERAGE_STEP <- 0.01' in main
assert 'monthly_plot <- monthly_plot_raw' in main
assert 'figure8' not in supp.lower() and 'figure9' not in supp.lower(), 'old plot names detected'
assert 'figure_ABC <- ggplot(' in province and province.count('ggsave(')==2
assert 'D:/Exel' not in province
assert all(f'"Figure_{f}_' in figs for f in (1,3,4,5))
with (root/'docs/figure_mapping.csv').open() as fh:
    rows=list(csv.DictReader(fh))
assert set(row['manuscript_figure'] for row in rows)=={'1','2','3','4','5','S1'}
print('PASS: mapping and structural checks (R not executed).')


# Check the author-supplied input and the two derived tables exist.
assert (root / 'data/bednet_distribution_by_province_2008_2022.xlsx').is_file()
with (root / 'data/province_campaigns_2009_2024.csv').open(newline='', encoding='utf-8') as fh:
    longrows=list(csv.DictReader(fh))
assert len(longrows)==26*16
assert len(set(r['province'] for r in longrows))==26
assert sorted(set(int(r['year']) for r in longrows))==list(range(2009,2025))
assert all(r['distributed'] in ('0','1') for r in longrows)
assert sum(int(r['distributed']) for r in longrows)==103
with (root / 'data/campaign_events_2009_2024.csv').open(newline='', encoding='utf-8') as fh:
    events=list(csv.DictReader(fh))
assert len(events)==103
print('PASS: Figure 2 data included; 26 provinces, 2009–2024, 103 event flags.')
