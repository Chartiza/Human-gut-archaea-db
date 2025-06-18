import os
import pandas as pd

# Directory path
dir_path = '/abud_results_arch21/'

# Get a list of all files in the directory that start with 'abud'
filenames = [f for f in os.listdir(dir_path) if f.startswith('abud')]

# Initialize an empty DataFrame to store the merged data
merged_df = pd.DataFrame()

for i, filename in enumerate(filenames):
    file_path = os.path.join(dir_path, filename)
    
    # Check if the file exists and is not empty
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        try:
            # Read each file into a pandas DataFrame
            df = pd.read_csv(file_path, sep='\t')
            
            # Ensure the DataFrame is not empty and has the expected structure
            if df.empty or df.shape[1] != 2:
                print(f"Skipping file {filename}: Empty or unexpected structure")
                continue
            
            # Get the file name without extension
            filename_no_ext = os.path.splitext(filename)[0]

            # If this is the first file, initialize the merged_df DataFrame
            if i == 0:
                merged_df = df
                merged_df.columns = ['rep_MAG_ID', filename_no_ext]  # Rename the column to the file's name (without extension)
            else:
                # If this is not the first file, merge it with the existing DataFrame
                df.columns = ['rep_MAG_ID', filename_no_ext]  # Rename the column to the file's name (without extension)
                merged_df = pd.merge(merged_df, df, on='rep_MAG_ID', how='outer')
        except pd.errors.EmptyDataError:
            print(f"Skipping file {filename}: Empty data")
        except pd.errors.ParserError:
            print(f"Skipping file {filename}: Parsing error")
        except Exception as e:
            print(f"Skipping file {filename}: Unexpected error {e}")
    else:
        print(f"Skipping file {filename}: Does not exist or is empty")

# Replace all NaN values with 0 (you can skip this step if NaN values are fine)
merged_df.fillna(0, inplace=True)

# Write the merged DataFrame to a CSV file
output_path = 'CoverM_merged_output_Arch21-Ill.csv'
merged_df.to_csv(output_path, index=False)

print(f"Merged data saved to {output_path}")
