#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第四章OCR识别脚本 - 使用按章节拆分的PDF
"""
import os
from pix2text import Pix2Text

pdf_path = "01_数学物理方程_按章节拆分/第04章_分离变量法_PDF页79-134.pdf"
output_dir = "01_数学物理方程_OCR_MD"

os.makedirs(output_dir, exist_ok=True)

print("Initializing Pix2Text...")
p2t = Pix2Text.from_config(device='cpu')

print(f"Starting OCR for {pdf_path}...")
try:
    result = p2t.recognize_pdf(
        pdf_path,
        auto_line_break=True
    )

    output_file = os.path.join(output_dir, "第04章_分离变量法_重新识别.md")
    print(f"Saving result to {output_file}...")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(result)

    print("OCR completed successfully!")
    print(f"Output saved to: {output_file}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
