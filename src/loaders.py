import cooler
import os
from typing import Optional

def get_cooler(filename: str, resolution: int = 1000000) -> cooler.Cooler:
    """
    Loads a cooler object at a specific resolution from the data/raw directory.
    
    Args:
        filename (str): Name of the file in data/raw (e.g., 'test.mcool')
        resolution (int): Resolution in bp (default 1MB). 
                          Common options: 1000000, 100000, 10000.
    
    Returns:
        cooler.Cooler: The loaded cooler object.
    
    Raises:
        FileNotFoundError: If the file does not exist in data/raw.
    """
    # Construct absolute path relative to this file
    # This assumes src/loaders.py structure, so we go up one level to root
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    file_path = os.path.join(base_dir, 'data', 'raw', filename)
    
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}. Did you run the download script?")
    
    # Check if it's an mcool (multi-resolution) file
    if filename.endswith('.mcool'):
        uri = f'{file_path}::resolutions/{resolution}'
    else:
        # Assume .cool files are single resolution and don't need uri suffix
        uri = file_path
    
    return cooler.Cooler(uri)

# Alias for backward compatibility if needed, or we can just use get_cooler
load_cooler = get_cooler
