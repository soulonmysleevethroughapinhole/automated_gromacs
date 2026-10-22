import os
import sys
import glob
import tarfile
import shutil
import pymol
from pymol import cmd

def main():
    if len(sys.argv) != 3:
        print("Usage: python render_pymol_movie.py <input_tar_gz> <output_mp4>", file=sys.stderr)
        sys.exit(1)

    tar_path = sys.argv[1]
    output_movie = sys.argv[2]
    
    # Dynamically find the working directory based on the output file path
    extract_dir = os.path.dirname(output_movie)
    frames_dir = os.path.join(extract_dir, "FRAMES")

    print(f"=== Starting Standalone PyMOL Renderer ===")
    print(f"[+] Extracting compressed coordinates from: {tar_path}")
    
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=extract_dir)

    print("[+] Initializing headless PyMOL...")
    pymol.finish_launching(['pymol', '-cq'])

    frame_files = sorted(glob.glob(os.path.join(frames_dir, "frame*.pdb")))
    print(f"[+] Found {len(frame_files)} frames. Loading into memory...")
    
    for f in frame_files:
        cmd.load(f, "trajectory_ensemble")

    cmd.show_as("cartoon", "all")
    cmd.color("cyan", "all")
    cmd.orient()

    print(f"[->] Encoding frames to: {output_movie}")
    cmd.movie.produce(output_movie, quality=90, width=1280, height=720)
    cmd.quit()

    print("[+] Cleaning up temporary frame files...")
    if os.path.exists(frames_dir):
        shutil.rmtree(frames_dir)
        
    print("=== PyMOL Movie Generation Complete ===")

if __name__ == "__main__":
    main()