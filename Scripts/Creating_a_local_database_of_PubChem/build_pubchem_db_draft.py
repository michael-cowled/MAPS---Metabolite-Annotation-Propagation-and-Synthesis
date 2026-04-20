import pandas as pd
import sqlite3
import os
import sys
import psutil

# --- Configuration ---
# Define the expected column names for each file type
COLUMN_MAP = {
    "CID-IUPAC": ["CID", "IUPAC"],
    "CID-Mass": ["CID", "Formula", "Monoisotopic.Mass", "Exact.Mass"],
    "CID-Title": ["CID", "Title"],
    "CID-SMILES": ["CID", "SMILES"],
    "CID-HMDB": ["CID", "Primary.HMDB.ID", "Secondary.HMDB.ID"] 
}

# List of files that DO NOT have a header row
NO_HEADER_FILES = ["CID-IUPAC", "CID-Mass", "CID-Title", "CID-SMILES"]

# Pause logic: How many chunks to process before asking to continue
CHUNKS_PER_BATCH = 10 
# ---------------------

DATA_DIR = os.getcwd()
DB_PATH = os.path.join(DATA_DIR, "PubChem_Indexed.sqlite")

def calculate_safe_chunksize():
    """Calculates a safe pandas chunksize based on currently available RAM."""
    available_mb = psutil.virtual_memory().available / (1024 * 1024)
    safe_mb = available_mb * 0.10  # Use only 10% of available free RAM
    estimated_chunksize = int(safe_mb * 1000) # ~1000 rows per MB estimate
    return max(10_000, min(estimated_chunksize, 500_000))

print(f"📁 Processing files in: {DATA_DIR}")
print(f"🗄️ Output database: {DB_PATH}")
print("-" * 30)

# 1. Discover all files
file_paths = {}
try:
    for filename in os.listdir(DATA_DIR):
        if os.path.isdir(os.path.join(DATA_DIR, filename)) or filename.endswith('.py') or filename.endswith('.sqlite') or filename.startswith('.'):
            continue
            
        path = os.path.join(DATA_DIR, filename)
        name = os.path.splitext(filename)[0]
        
        if name in COLUMN_MAP:
            file_paths[name] = path
        else:
            print(f"   ⚠️ Skipping {filename}: No column map defined.")

except Exception as e:
    print(f"\n❌ Error listing files: {e}")
    sys.exit(1)

if not file_paths:
    print("🛑 No valid files were found. Exiting.")
    sys.exit(1)

# 2. Determine optimal RAM chunksize
chunksize = calculate_safe_chunksize()
print(f"\n🧠 Available RAM optimized: Processing files in chunks of {chunksize:,} rows.")

# 3. Load files into separate SQLite tables (with Resume Logic)
conn = sqlite3.connect(DB_PATH)
print(f"\n📄 Found {len(file_paths)} files. Loading into SQLite...")

for name, path in file_paths.items():
    table_name = name.replace('-', '_')
    print(f"\n   -> Processing {name} into table '{table_name}'...")
    
    col_names = COLUMN_MAP[name]
    has_header = name not in NO_HEADER_FILES
    
    # --- Check Database for Existing Progress ---
    rows_already_processed = 0
    cursor = conn.execute(f"SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{table_name}'")
    if cursor.fetchone()[0] == 1:
        cursor = conn.execute(f"SELECT COUNT(*) FROM {table_name}")
        rows_already_processed = cursor.fetchone()[0]
        if rows_already_processed > 0:
            print(f"      Resume detected! Skipping the first {rows_already_processed:,} rows.")
    
    # --- Calculate Skiprows ---
    if has_header and rows_already_processed > 0:
        # Keep the header (line 0), skip the data lines already processed
        skip_logic = range(1, rows_already_processed + 1) 
        read_kwargs = {"header": 0}
    else:
        # No header, just skip the exact number of rows we've already done
        skip_logic = rows_already_processed
        read_kwargs = {"header": None, "names": col_names}
    
    # --- Process the File in Pausable Batches ---
    try:
        chunk_iterator = pd.read_csv(
            path, 
            sep="\t", 
            dtype={"CID": str},
            chunksize=chunksize,
            skiprows=skip_logic,
            **read_kwargs
        )
        
        chunks_processed_this_session = 0
        
        for i, chunk in enumerate(chunk_iterator):
            if 'CID' not in chunk.columns:
                 print(f"      ❌ Skipping {name}: 'CID' column missing.")
                 break
            
            write_mode = "replace" if (rows_already_processed == 0 and chunks_processed_this_session == 0) else "append"
            
            chunk.to_sql(table_name, conn, if_exists=write_mode, index=False)
            chunks_processed_this_session += 1
            
            total_current_rows = rows_already_processed + (chunks_processed_this_session * chunksize)
            print(f"      ... Processed roughly {total_current_rows:,} rows total")

            # Pause Mechanism
            if chunks_processed_this_session % CHUNKS_PER_BATCH == 0:
                user_choice = input(f"\n⏸️ Paused after {CHUNKS_PER_BATCH} chunks. Continue? (y/n): ")
                if user_choice.lower() != 'y':
                    print("🛑 Stopping execution. Your progress is saved in the database.")
                    conn.close()
                    sys.exit(0)

        # Create index only after the entire file is 100% finished
        print(f"      ✅ Finished loading {name}. Creating index...")
        conn.execute(f"CREATE INDEX IF NOT EXISTS idx_{table_name}_cid ON {table_name} (CID);")
            
    except pd.errors.EmptyDataError:
        print(f"      ✅ File {name} is already completely processed.")
    except Exception as e:
        print(f"      ❌ Failed to load {name}: {e}")

# 4. Perform the massive merge on DISK using SQLite
# Note: This single operation cannot be easily paused. Let it run to completion.
print("\n🤝 Finalizing: Merging data on disk via SQLite...")
print("   (This step might take a while and should not be interrupted. Ensure sufficient disk space.)")

# Step 4a: Create a master list of all unique CIDs
print("   -> Compiling master list of unique CIDs...")
union_queries = [f"SELECT CID FROM {name.replace('-', '_')}" for name in file_paths.keys()]
union_sql = "\nUNION\n".join(union_queries)
conn.execute("DROP TABLE IF EXISTS master_cids;")
conn.execute(f"CREATE TABLE master_cids AS {union_sql};")
conn.execute("CREATE INDEX idx_master_cid ON master_cids (CID);")

# Step 4b: Construct and execute the LEFT JOIN query
print("   -> Executing final join into 'pubchem_data'...")
select_cols = ["m.CID"]
joins = []

for name in file_paths.keys():
    table_name = name.replace('-', '_')
    cols = [c for c in COLUMN_MAP[name] if c != "CID"]
    for c in cols:
        select_cols.append(f"{table_name}.\"{c}\"") 
    
    joins.append(f"LEFT JOIN {table_name} ON m.CID = {table_name}.CID")

select_sql = ",\n        ".join(select_cols)
join_sql = "\n    ".join(joins)

final_join_query = f"""
CREATE TABLE pubchem_data AS
SELECT 
    {select_sql}
FROM master_cids m
    {join_sql};
"""
conn.execute("DROP TABLE IF EXISTS pubchem_data;")
conn.execute(final_join_query)

# 5. Final Index and Cleanup
print("🔑 Indexing the final combined table...")
conn.execute("CREATE INDEX idx_final_cid ON pubchem_data (CID);")

print("🧹 Cleaning up temporary tables to save hard drive space...")
conn.execute("DROP TABLE master_cids;")
for name in file_paths.keys():
    table_name = name.replace('-', '_')
    conn.execute(f"DROP TABLE {table_name};")

print("🗜️ Vacuuming database (defragmenting and reclaiming space)...")
conn.execute("VACUUM;")
conn.close()

print(f"\n✅ Success! Database safely built at: {DB_PATH}")