#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第三章OCR识别脚本
"""
import os
from pix2text import Pix2Text

pdf_path = "/Users/tmo/Documents/数学物理方程/课本教材和学习指导/01_数学物理方程_按章节拆分/第03章_波动方程的初值柯西问题与行波法_PDF页52-78.pdf"
output_dir = "/Users/tmo/Documents/数学物理方程/课本教材和学习指导/01_数学物理方程_OCR_MD"

os.makedirs(output_dir, exist_ok=True)

print("Initializing Pix2Text...")
# 关键：使用from_config并指定device='cpu'
p2t = Pix2Text.from_config(device='cpu')

print(f"Starting OCR for {pdf_path}...")
try:
    result = p2t.recognize_pdf(
        pdf_path,
        auto_line_break=True
    )

    output_file = os.path.join(output_dir, "第03章_波动方程的初值柯西问题与行波法.md")
    print(f"Saving result to {output_file}...")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(result)

    print("OCR completed successfully!")
    print(f"Output saved to: {output_file}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
