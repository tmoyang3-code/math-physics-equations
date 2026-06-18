#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys

# 设置环境变量
os.environ['OCR_PROVIDER'] = 'cpu'

from pix2text import Pix2Text

# PDF文件路径
pdf_path = "/Users/tmo/Documents/数学物理方程/课本教材和学习指导/01_数学物理方程_按章节拆分/第03章_波动方程的初值柯西问题与行波法_PDF页52-78.pdf"
output_dir = "/Users/tmo/Documents/数学物理方程/课本教材和学习指导/01_数学物理方程_OCR_MD"
output_md = os.path.join(output_dir, "第03章_波动方程的初值柯西问题与行波法.md")

# 确保输出目录存在
os.makedirs(output_dir, exist_ok=True)

print("Initializing Pix2Text...")
p2t = Pix2Text.from_config(device='cpu')

print(f"Starting OCR for PDF: {pdf_path}")
print("This may take 15-30 minutes for a 27-page document...")

try:
    # 识别整个PDF
    result = p2t.recognize_pdf(pdf_path, auto_line_break=True)

    print("Recognition completed!")
    print("Result type:", type(result))

    # 尝试不同的方式获取内容
    md_content = ""

    # 方式1: 尝试直接转换为字符串
    try:
        md_content = str(result)
        print("Method 1: str() worked, length:", len(md_content))
    except Exception as e:
        print(f"Method 1 failed: {e}")

    # 方式2: 尝试 to_markdown() 方法
    if len(md_content) < 100 and hasattr(result, 'to_markdown'):
        try:
            md_content = result.to_markdown()
            print("Method 2: to_markdown() worked, length:", len(md_content))
        except Exception as e:
            print(f"Method 2 failed: {e}")

    # 方式3: 尝试 pages 属性
    if len(md_content) < 100 and hasattr(result, 'pages'):
        try:
            print(f"Found {len(result.pages)} pages")
            for i, page in enumerate(result.pages):
                print(f"Processing page {i}...")
                if hasattr(page, 'to_markdown'):
                    page_md = page.to_markdown()
                    md_content += f"\n\n---\n## 第 {i+1} 页\n\n{page_md}\n"
                else:
                    md_content += f"\n\n---\n## 第 {i+1} 页\n\n{str(page)}\n"
            print("Method 3: pages to_markdown worked, total length:", len(md_content))
        except Exception as e:
            print(f"Method 3 failed: {e}")

    print(f"\n\nTotal markdown content length: {len(md_content)}")

    # 保存结果
    with open(output_md, 'w', encoding='utf-8') as f:
        f.write(md_content)

    print(f"\n✓ Markdown file saved to: {output_md}")

    # 预览前几行
    print("\nPreview (first 500 characters):")
    print(md_content[:500])

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
