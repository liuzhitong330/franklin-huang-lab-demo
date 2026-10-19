"""Wrap the checked JSON for the static page. Does not change any values."""
import argparse
import json
from pathlib import Path

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    data=json.loads(args.data.read_text())
    text=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    args.output.write_text("(function () { 'use strict'; window.HUANG_DATA = "+text+"; })();\n")
