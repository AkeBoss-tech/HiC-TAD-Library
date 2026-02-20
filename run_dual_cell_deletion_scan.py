#!/usr/bin/env python3
"""
Run deletion scan for both mouse cell types.

Usage:
    python run_dual_cell_deletion_scan.py jingyun --del-sizes 10 40 80
    python run_dual_cell_deletion_scan.py edward --del-sizes 10 40 80
"""
import os
import sys
import subprocess

# Cell types to test
CELL_TYPES = [
    ('EFO:0004038', 'Mouse_ESC'),
    ('CL:0000207', 'Olfactory_receptor_cell'),
]

def run_scan_for_cell_type(region, del_sizes, cell_type_id, cell_name):
    """Run deletion_scan.py with a specific cell type."""
    print(f'\n{"="*70}')
    print(f'Running deletion scan for: {cell_name} ({cell_type_id})')
    print(f'{"="*70}\n')

    # Read the original script
    with open('deletion_scan.py', 'r') as f:
        original_content = f.read()

    # Modify cell type
    modified_content = original_content.replace(
        "CELL_TYPE = 'EFO:0004038'",
        f"CELL_TYPE = '{cell_type_id}'"
    ).replace(
        "CELL_NAME = 'Mouse ESC'",
        f"CELL_NAME = '{cell_name.replace('_', ' ')}'"
    ).replace(
        "SAFE         = cfg['safe']",
        f"SAFE         = cfg['safe'] + '_{cell_name.lower()}'"
    )

    # Write temporary modified script
    temp_script = f'deletion_scan_temp_{cell_name}.py'
    with open(temp_script, 'w') as f:
        f.write(modified_content)

    # Run it
    try:
        cmd = ['python', temp_script, region]
        if del_sizes:
            cmd.extend(['--del-sizes'] + [str(s) for s in del_sizes])

        subprocess.run(cmd, check=True)

        print(f'\n✓ Completed scan for {cell_name}')

    finally:
        # Clean up temp script
        if os.path.exists(temp_script):
            os.remove(temp_script)


def main():
    if len(sys.argv) < 2:
        print("Usage: python run_dual_cell_deletion_scan.py <region> [--del-sizes SIZE1 SIZE2 ...]")
        print("Example: python run_dual_cell_deletion_scan.py jingyun --del-sizes 10 40 80")
        sys.exit(1)

    region = sys.argv[1]

    # Parse deletion sizes
    del_sizes = None
    if '--del-sizes' in sys.argv:
        idx = sys.argv.index('--del-sizes')
        del_sizes = [int(x) for x in sys.argv[idx+1:]]

    print("="*70)
    print("Multi-Cell Type Deletion Scan")
    print("="*70)
    print(f"Region: {region}")
    print(f"Deletion sizes: {del_sizes if del_sizes else 'original insulator size'}")
    print(f"Cell types: {len(CELL_TYPES)}")
    print("="*70)

    # Run for each cell type
    for cell_type_id, cell_name in CELL_TYPES:
        run_scan_for_cell_type(region, del_sizes, cell_type_id, cell_name)

    print("\n" + "="*70)
    print("ALL CELL TYPES COMPLETE!")
    print("="*70)
    print("Check media/ directory for all visualizations")
    print("HTML reports generated for each cell type")


if __name__ == '__main__':
    main()
