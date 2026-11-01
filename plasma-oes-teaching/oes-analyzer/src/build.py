"""Rebuild oes-analyzer/index.html from the sources in this folder.

Usage (from this folder):
    npm pack plotly.js-dist-min@2.27.0 && tar xzf plotly.js-dist-min-2.27.0.tgz
    python3 build.py

Plotly is inlined so the finished page works offline.
"""
from pathlib import Path
here = Path(__file__).parent
ui = (here / 'ui.html').read_text(encoding='utf-8')
eng = (here / 'engine.js').read_text(encoding='utf-8')
app = (here / 'app.js').read_text(encoding='utf-8')
pl = (here / 'package' / 'plotly.min.js').read_text(encoding='utf-8')
pl = pl.replace('</script', '<\\/script').replace('�', '\\ufffd')
head = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">'
        '<style>:root{color-scheme:light}body{margin:0;font:14px system-ui,sans-serif}'
        'img{max-width:100%}[hidden]{display:none!important}</style></head><body>\n')
page = head + ui + '\n<script>' + pl + '</script>\n<script>\n' + eng + '\n</script>\n<script>\n' + app + '\n</script>\n</body></html>\n'
(here.parent / 'index.html').write_text(page, encoding='utf-8')
print('wrote', here.parent / 'index.html', f'({len(page)/1e6:.1f} MB)')
