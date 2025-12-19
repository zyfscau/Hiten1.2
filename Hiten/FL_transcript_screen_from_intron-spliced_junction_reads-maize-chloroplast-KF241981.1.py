import pandas as pd

# Read CSV file
df = pd.read_csv("vlookup-result-sift.csv")

# Define 5'end ranges
five_prime_cond = df.iloc[:, 73].between(69471, 69483) | df.iloc[:, 73].between(70313, 70331)

# Define 3'end ranges (combine all possible intervals)
three_prime_cond = (
    df.iloc[:, 74].between(69209, 69221) |
    df.iloc[:, 74].between(68696, 68711) |
    df.iloc[:, 74].between(91474, 91494)
)

# Define strand condition
strand_cond = (df.iloc[:, 76] == '-')

# Combine all conditions
filtered_df = df[five_prime_cond & three_prime_cond & strand_cond]

# Save results to new CSV file
filtered_df.to_csv("intron-FL.csv", index=False)

# Output statistics
print(f"Filtering completed! Number of records matching criteria: {len(filtered_df)}")
print(f"5'end ranges: 69471-69483 or 70313-70331")
print(f"3'end ranges: 69209-69221, 68696-68711, 91474-91494")
print(f"Strand: Negative strand (-)")
print(f"Results saved to: intron-FL.csv")