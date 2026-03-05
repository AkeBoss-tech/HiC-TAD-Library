import os
import markdown
import glob

def build_docs():
    """Convert Markdown documentation files to HTML."""
    docs_files = [
        "README.md",
        "ANALYSIS.md",
        "TESTING.md",
        "docs/notebooks.md"
    ]

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title} - HiC-TAD Library</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 900px;
            margin: 0 auto;
            padding: 40px 20px;
            background-color: #f9f9f9;
        }}
        .container {{
            background: white;
            padding: 40px;
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        h1, h2, h3, h4, h5, h6 {{
            color: #2196F3;
            margin-top: 1.5em;
            margin-bottom: 0.5em;
        }}
        h1 {{ border-bottom: 2px solid #eee; padding-bottom: 10px; }}
        h2 {{ border-bottom: 1px solid #eee; padding-bottom: 5px; }}
        a {{ color: #1976D2; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 5px;
            border-radius: 3px;
            font-family: Consolas, Monaco, 'Courier New', monospace;
            font-size: 0.9em;
        }}
        pre {{
            background-color: #282c34;
            color: #abb2bf;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            line-height: 1.4;
        }}
        pre code {{
            background-color: transparent;
            padding: 0;
            color: inherit;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 10px;
            text-align: left;
        }}
        th {{
            background-color: #f2f2f2;
            color: #333;
        }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        blockquote {{
            border-left: 4px solid #2196F3;
            margin: 0;
            padding-left: 15px;
            color: #666;
            background-color: #f5fafe;
            padding: 10px 15px;
            border-radius: 0 4px 4px 0;
        }}
        img {{
            max-width: 100%;
            height: auto;
            display: block;
            margin: 20px auto;
            border-radius: 4px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .back-link {{
            display: inline-block;
            margin-bottom: 20px;
            margin-top: -10px;
            padding: 8px 15px;
            background-color: #f0f0f0;
            border-radius: 4px;
            font-weight: bold;
            transition: background-color 0.2s;
        }}
        .back-link:hover {{ background-color: #e0e0e0; text-decoration: none; }}
    </style>
</head>
<body>
    <div class="container">
        <a href="index.html" class="back-link">← Back to Gallery</a>
        {content}
    </div>
</body>
</html>"""

    # Ensure output directory exists if needed, but we'll save to root for easy linking
    
    for md_file in docs_files:
        if not os.path.exists(md_file):
            print(f"Warning: {md_file} not found. Skipping.")
            continue
            
        with open(md_file, "r", encoding="utf-8") as f:
            text = f.read()
            
        # Convert Markdown to HTML
        # Using extensions for tables, fenced code blocks, and TOC
        html_content = markdown.markdown(
            text, 
            extensions=['fenced_code', 'tables', 'toc', 'sane_lists']
        )
        
        # Get filename without path and extension
        base_name = os.path.basename(md_file)
        name_only = os.path.splitext(base_name)[0]
        
        # Output file path (always in root for easy linking from index.html)
        out_file = f"{name_only}.html"
        
        # Create full HTML page
        full_html = html_template.format(
            title=name_only.replace("_", " ").title(),
            content=html_content
        )
        
        with open(out_file, "w", encoding="utf-8") as f:
            f.write(full_html)
            
        print(f"Built {out_file} from {md_file}")

if __name__ == "__main__":
    build_docs()
