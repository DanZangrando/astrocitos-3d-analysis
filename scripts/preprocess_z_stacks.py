
import numpy as np
import tifffile
from pathlib import Path
from tqdm import tqdm
import shutil
import argparse
import os

def create_lut(color):
    """
    Creates a 768-element numpy array representing an ImageJ LUT (RGB).
    """
    r = np.zeros(256, dtype=np.uint8)
    g = np.zeros(256, dtype=np.uint8)
    b = np.zeros(256, dtype=np.uint8)
    
    ramp = np.arange(256, dtype=np.uint8)
    
    if color == 'blue':
        b = ramp
    elif color == 'green':
        g = ramp
    elif color == 'red':
        r = ramp
    elif color == 'cyan':
        g = ramp
        b = ramp
    elif color == 'magenta':
        r = ramp
        b = ramp
    elif color == 'yellow':
        r = ramp
        g = ramp
    elif color == 'gray':
        r = ramp
        g = ramp
        b = ramp
        
    return np.concatenate([r, g, b]).astype(np.uint8)

def get_best_10_z_slices(image, gfap_channel_index=1):
    """
    Identifies the start index of the 10 continuous Z-slices with the maximum sum intensity
    in the GFAP channel.
    """
    if image.ndim != 4:
        raise ValueError(f"Expected 4D image, got {image.ndim}D")
    
    # Ensure (Z, C, Y, X)
    if image.shape[0] < image.shape[1]:
        image = image.transpose(1, 0, 2, 3)
        
    n_z = image.shape[0]
    n_c = image.shape[1]
    
    if n_z < 10:
        return -1, image
        
    if gfap_channel_index >= n_c:
        gfap_channel_index = 1 if n_c > 1 else 0
            
    gfap_data = image[:, gfap_channel_index, :, :]
    z_intensities = [np.sum(gfap_data[i]) for i in range(n_z)]
    
    window_size = 10
    max_sum = -1.0
    best_start = 0
    
    for i in range(n_z - window_size + 1):
        current_sum = sum(z_intensities[i : i + window_size])
        if current_sum > max_sum:
            max_sum = current_sum
            best_start = i
            
    return best_start, image

def save_image_formatted(output_path, image):
    """
    Saves the image with only 2 channels (DAPI, GFAP) and correct colors.
    """
    if image.ndim != 4:
        return False
        
    # Trim to 2 channels if more than 2
    if image.shape[1] > 2:
        image = image[:, :2, :, :]
    
    if image.shape[1] < 1:
        return False
        
    # Prepare Metadata
    ij_metadata = {
        'axes': 'ZCYX',
        'LUTs': [] 
    }
    
    luts = []
    if image.shape[1] >= 1:
        luts.append(create_lut('blue')) # Ch0
    if image.shape[1] >= 2:
        luts.append(create_lut('green')) # Ch1
        
    ij_metadata['LUTs'] = luts
    
    # Use TiffWriter for better control
    with tifffile.TiffWriter(str(output_path), imagej=True) as tif:
        tif.write(image, metadata=ij_metadata)
        
    return True

def process_file(file_path, output_root, input_root):
    try:
        rel_path = file_path.relative_to(input_root)
        output_path = output_root / rel_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with tifffile.TiffFile(str(file_path)) as tf:
            image = tf.asarray()
            
        if image.ndim == 3:
             image = np.expand_dims(image, axis=1) # (Z, 1, Y, X)
             
        if image.ndim == 4:
            s = image.shape
            if s[0] < s[1]:
                 image = np.transpose(image, (1, 0, 2, 3))
        
        n_z = image.shape[0]
        
        if n_z < 10:
            print(f"WARNING: {file_path.name} has fewer than 10 Z-slices ({n_z}). Saving as is (channels trimmed).")
            # Always format correctly
            save_image_formatted(output_path, image)
            return "skipped_short"
            
        else:
            if n_z == 10:
                cropped_image = image
                result_type = "fixed_colors"
            else:
                start_idx, _ = get_best_10_z_slices(image)
                cropped_image = image[start_idx : start_idx + 10, ...]
                result_type = "processed"
            
            save_image_formatted(output_path, cropped_image)
            return result_type

    except Exception as e:
        print(f"ERROR processing {file_path}: {e}")
        return "error"

def main():
    parser = argparse.ArgumentParser(description="Filter Z-stacks to best 10 continuous slices and fix colors.")
    parser.add_argument("--input", default="data/raw", help="Input directory (default: data/raw)")
    parser.add_argument("--output", default="data/raw_z_corrected", help="Output directory (default: data/raw_z_corrected)")
    
    args = parser.parse_args()
    
    input_root = Path(args.input)
    output_root = Path(args.output)
    
    if not input_root.exists():
        print(f"Input directory not found: {input_root}")
        return

    files = sorted(list(input_root.rglob("*.tif")) + list(input_root.rglob("*.tiff")))
    
    print(f"Found {len(files)} images in {input_root}")
    print(f"Saving processed images to {output_root}")
    
    stats = {"processed": 0, "fixed_colors": 0, "skipped_short": 0, "error": 0}
    
    for p in tqdm(files):
        parts = [part.lower() for part in p.parts]
        res = process_file(p, output_root, input_root)
        stats[res] += 1
        
    print("\n--- Processing Summary ---")
    print(f"Total files found: {len(files)}")
    print(f"Processed (cropped < 10z): {stats['processed']}")
    print(f"Fixed Colors (already 10z): {stats['fixed_colors']}")
    print(f"Skipped Short (< 10z): {stats['skipped_short']}")
    print(f"Errors: {stats['error']}")
    print(f"\nAll operations complete. Output is in {output_root}")

if __name__ == "__main__":
    main()
