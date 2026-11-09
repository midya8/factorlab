import argparse,json,sys
from pathlib import Path
from .core import audit,strict_loads

def main(argv=None):
    p=argparse.ArgumentParser(description='Exact binary factor elimination and mixed sum/max audits');p.add_argument('input',type=Path);p.add_argument('--output',type=Path);a=p.parse_args(argv)
    try:r=audit(strict_loads(a.input.read_text()));text=json.dumps(r,indent=2)+'\n'
    except (ValueError,OSError,TypeError) as e:print(str(e),file=sys.stderr);return 2
    if a.output:a.output.write_text(text)
    else:print(text,end='')
    return 0
