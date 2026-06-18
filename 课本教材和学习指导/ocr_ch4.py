#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第四章OCR识别脚本
"""
import os
from pix2text import Pix2Text

pdf_path = "01_数学物理方程 (陈才生) (z-library.sk, 1lib.sk, z-lib.sk).pdf"
output_dir = "01_数学物理方程_OCR_MD/ch4_revised"

os.makedirs(output_dir, exist_ok=True)

print("Initializing Pix2Text...")
p2t = Pix2Text.from_config(device='cpu')

print(f"Starting OCR for pages 67-118...")
try:
    result = p2t.recognize_pdf(
        pdf_path,
        start_page=67,
        end_page=118,
        auto_line_break=True
    )

    output_file = os.path.join(output_dir, "第04章_分离变量法_重新识别.md")
    print(f"Saving result to {output_file}...")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(result)

    print("OCR completed successfully!")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
