import pandas as pd

# Read CSV file
df = pd.read_csv("vlookup-result-sift.csv")

# Define 5'end ranges
five_prime_cond = df.iloc[:, 73].between(67562, 67574) | df.iloc[:, 73].between(68402, 68419)

# Define 3'end ranges (combine all possible intervals)
three_prime_cond = (
    df.iloc[:, 74].between(85255, 85269) |
    df.iloc[:, 74].between(85347, 85360) |
    df.iloc[:, 74].between(86964, 86978) |
    df.iloc[:, 74].between(87641, 87654) |
    df.iloc[:, 74].between(87788, 87820)
)

# Define strand condition
strand_cond = (df.iloc[:, 76] == '-')

# Combine all conditions
filtered_df = df[five_prime_cond & three_prime_cond & strand_cond]

# Save results to new CSV file
filtered_df.to_csv("intron-FL.csv", index=False)

# Output statistics
print(f"Filtering completed! Number of records matching criteria: {len(filtered_df)}")
print(f"5'end ranges: 67562-67574 or 68402-68419")
print(f"3'end ranges: 85255-85269, 85347-85360, 86964-86978, 87641-87654, 87788-87820")
print(f"Strand: Negative strand (-)")
print(f"Results saved to: intron-FL.csv")