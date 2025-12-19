import pandas as pd
import sqlite3
import os
import sys
from functools import reduce

# --- Configuration ---
# Define the expected column names for each file type
COLUMN_MAP = {
    "CID-IUPAC": ["CID", "IUPAC"],
    "CID-Mass": ["CID", "Formula", "Monoisotopic.Mass", "Exact.Mass"],
    "CID-Title": ["CID", "Title"],
    "CID-SMILES": ["CID", "SMILES"],
    "CID-HMDB": ["CID", "Primary.HMDB.ID", "Secondary.HMDB.ID"] 
}

# List of files that still DO NOT have a header row
NO_HEADER_FILES = ["CID-IUPAC", "CID-Mass", "CID-Title", "CID-SMILES"]
# ---------------------

# 1. Set DATA_DIR to the current working directory
DATA_DIR = os.getcwd()
DB_PATH = os.path.join(DATA_DIR, "PubChem_Indexed.sqlite")

print(f"📁 Processing files in current working directory: {DATA_DIR}")
print(f"🗄️ Output database will be: {DB_PATH}")
print("-" * 30)

# 2. Discover all files in the directory
file_paths = {}
try:
    for filename in os.listdir(DATA_DIR):
        # Skip directories and non-data files
        if os.path.isdir(os.path.join(DATA_DIR, filename)) or filename.endswith('.py') or filename.endswith('.sqlite') or filename.startswith('.'):
            continue
            
        path = os.path.join(DATA_DIR, filename)
        name = os.path.splitext(filename)[0] # e.g., "CID-IUPAC"
        
        # Only process files whose names are in our COLUMN_MAP
        if name in COLUMN_MAP:
            file_paths[name] = path
        else:
            print(f"   ⚠️ Skipping {filename}: No column map defined for this file type.")

except Exception as e:
    print(f"\nAn error occurred while listing files: {e}")
    sys.exit(1)


# 3. Load all discovered files
dataframes = {}
print(f"\n📄 Found {len(file_paths)} files to process...")

for name, path in file_paths.items():
    print(f"   Loading and renaming columns for {name}...")
    try:
        col_names = COLUMN_MAP[name]
        
        # Determine if the file has a header or not
        if name in NO_HEADER_FILES:
            # Original files (no header)
            read_kwargs = {"header": None, "names": col_names}
        else:
            # Fixed file (header included)
            # The 'names' argument is omitted here, as pandas will use the header row
            # which matches our col_names, and we avoid assigning headers twice.
            read_kwargs = {"header": 0} 
        
        df = pd.read_csv(
            path, 
            sep="\t", 
            dtype={"CID": str},
            **read_kwargs # Pass the dynamically selected arguments
        )
        
        # NOTE: The column count check is complicated by the header vs. no-header scenario.
        # We rely on the files being correct. The main thing is to ensure 'CID' is present.
        if 'CID' not in df.columns:
             print(f"      ❌ Skipping {name}: 'CID' column not found after loading.")
             continue
            
        dataframes[name] = df.set_index("CID")
        
    except Exception as e:
        print(f"      ❌ Failed to load {name}: {e}")

# 4. Merge all dataframes (unchanged)
if not dataframes:
    print("🛑 No valid files were loaded. Exiting.")
    sys.exit(1)

print("\n🤝 Merging dataframes...")
merged_df = reduce(lambda left, right: left.join(right, how="outer"), dataframes.values())
merged_df.reset_index(inplace=True)

# 5 & 6. Save to SQLite and index (unchanged)
print("\n💾 Saving to SQLite...")
conn = sqlite3.connect(DB_PATH)
merged_df.to_sql("pubchem_data", conn, if_exists="replace", index=False)

print("🔑 Creating index on CID...")
with conn:
    conn.execute("CREATE INDEX IF NOT EXISTS idx_cid ON pubchem_data (CID);")

conn.close()

print(f"\n✅ SQLite database created and indexed at: {DB_PATH}")
