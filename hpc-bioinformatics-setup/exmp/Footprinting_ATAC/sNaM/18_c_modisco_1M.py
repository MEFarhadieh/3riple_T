#!/usr/bin/env python3
"""
USAGE:
  python3 19_standalone_html.py \
    --html   18_modisco_1M/report/motifs.html \
    --jaspar 18_modisco_1M/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt \
    --out    18_modisco_1M/modisco_standalone.html \
    --base   /mnt/archive/farhadie/tn5_bias/skin_Mphage/sNaM_sorted

"""

import argparse
import base64
import os
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser()
parser.add_argument('--html',   required=True, help='HTML report from modisco')
parser.add_argument('--jaspar', required=True, help='JASPAR meme for  ID→name mapping')
parser.add_argument('--out',    required=True, help='HTML output standalone')
parser.add_argument('--base',   default=None,  help='base path for resolve image paths ( HTML)')
args = parser.parse_args()

BASE_PATH = args.base if args.base else os.path.dirname(os.path.abspath(args.html))

print("خواندن JASPAR names...")
jaspar_names = {}
with open(args.jaspar) as f:
    for line in f:
        if line.startswith('MOTIF'):
            parts = line.strip().split()
            if len(parts) >= 3:
                jaspar_names[parts[1]] = parts[2]
print(f"JASPAR motif number: {len(jaspar_names)}")

def img_to_base64(src, base_path):
    path = src.replace('file://', '')
    if not os.path.isabs(path):
        path = os.path.join(base_path, path)
    try:
        with open(path, 'rb') as f:
            data = base64.b64encode(f.read()).decode()
        return f"data:image/png;base64,{data}"
    except:
        return src

# ─── parse HTML ───
print("پردازش HTML...")
with open(args.html) as f:
    soup = BeautifulSoup(f.read(), 'html.parser')

table = soup.find('table')
thead = table.find('thead')
tbody = table.find('tbody')

# header
header_row = thead.find('tr')
for col in ['match0_TF', 'match1_TF', 'match2_TF']:
    th = soup.new_tag('th')
    th.string = col
    header_row.append(th)


for row in tbody.find_all('tr'):
    cells = row.find_all('td')
    for match_idx in [4, 7, 10]:
        td = soup.new_tag('td')
        if match_idx < len(cells):
            ma_id = cells[match_idx].get_text(strip=True)
            td.string = jaspar_names.get(ma_id, '?')
        row.append(td)

# ─── embed images ───
print("Images embeding...")
for img in soup.find_all('img'):
    src = img.get('src', '')
    img['src'] = img_to_base64(src, BASE_PATH)
    img['width'] = '200'

# ─── CSS ───
if not soup.find('head'):
    head = soup.new_tag('head')
    soup.insert(0, head)
soup.find('head').append(soup.new_tag('meta', charset='utf-8'))

style = soup.new_tag('style')
style.string = """
body { font-family: Arial, sans-serif; font-size: 11px; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #ccc; padding: 4px; text-align: center; vertical-align: middle; }
th { background: #4472C4; color: white; position: sticky; top: 0; }
tr:nth-child(even) { background: #f9f9f9; }
img { max-width: 200px; }
"""
soup.find('head').append(style)

# ─── saving ───
with open(args.out, 'w') as f:
    f.write(str(soup))

size = os.path.getsize(args.out) / 1e6
print(f"saved: {args.out}  ({size:.1f} MB)")
