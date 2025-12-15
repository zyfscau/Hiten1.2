import pandas as pd

# Read CSV file
df = pd.read_csv("vlookup-result-sift.csv")

# Define column indices
FIVE_PRIME_COL = 73   # Column containing 5'end
THREE_PRIME_COL = 74  # Column containing 3'end
STRAND_COL = 76       # Column containing strand

# Condition group 1: four sub-conditions
cond1_1 = (
    (df.iloc[:, FIVE_PRIME_COL].between(49267, 49280) |
    df.iloc[:, FIVE_PRIME_COL].between(320889, 320991) |
    df.iloc[:, FIVE_PRIME_COL].between(547952, 547962)
) & (
    df.iloc[:, THREE_PRIME_COL].between(52415, 52427) |
    df.iloc[:, THREE_PRIME_COL].between(323749, 323759) |
    df.iloc[:, THREE_PRIME_COL].between(549447, 549459)
) & (df.iloc[:, STRAND_COL] == '+')
) 
cond1_2 = (
    (df.iloc[:, FIVE_PRIME_COL].between(123741, 123753) |
    df.iloc[:, FIVE_PRIME_COL].between(267449, 267459))
) & (
    (df.iloc[:, THREE_PRIME_COL].between(266885, 266898) |
    df.iloc[:, THREE_PRIME_COL].between(120593, 120605))
) & (df.iloc[:, STRAND_COL] == '-')

cond1_3 = (
    (df.iloc[:, FIVE_PRIME_COL].between(49267, 49280) |
    df.iloc[:, FIVE_PRIME_COL].between(320889, 320991) |
    df.iloc[:, FIVE_PRIME_COL].between(547952, 547962))
) & (
    (df.iloc[:, THREE_PRIME_COL].between(266885, 266898) |
    df.iloc[:, THREE_PRIME_COL].between(120593, 120605))
) & (df.iloc[:, STRAND_COL] == '+/-')
cond1_4 = (
    (df.iloc[:, FIVE_PRIME_COL].between(123741, 123753) |
    df.iloc[:, FIVE_PRIME_COL].between(267449, 267459))
) & (
    (df.iloc[:, THREE_PRIME_COL].between(52415, 52427) |
    df.iloc[:, THREE_PRIME_COL].between(323749, 323759) |
    df.iloc[:, THREE_PRIME_COL].between(549447, 549459))
) & (df.iloc[:, STRAND_COL] == '-/+')

# Condition group 2: four sub-conditions
cond2_1 = (
    df.iloc[:, FIVE_PRIME_COL].between(111274, 111288)
) & (
    df.iloc[:, THREE_PRIME_COL].between(115856, 115868)
) & (df.iloc[:, STRAND_COL] == '+')

cond2_2 = (
    (df.iloc[:, FIVE_PRIME_COL].between(513364, 513375) |
    df.iloc[:, FIVE_PRIME_COL].between(512936, 512947))  # Corrected to ascending order
) & (
    (df.iloc[:, THREE_PRIME_COL].between(57152, 57164) |
    df.iloc[:, THREE_PRIME_COL].between(509329, 509340))
) & (df.iloc[:, STRAND_COL] == '-')

cond2_3 = (
    df.iloc[:, FIVE_PRIME_COL].between(111274, 111288)
) & (
    (df.iloc[:, THREE_PRIME_COL].between(57152, 57164) |
    df.iloc[:, THREE_PRIME_COL].between(509329, 509340))
) & (df.iloc[:, STRAND_COL] == '+/-')

cond2_4 = (
    (df.iloc[:, FIVE_PRIME_COL].between(513364, 513375) |
    df.iloc[:, FIVE_PRIME_COL].between(512936, 512947))  # Corrected to ascending order
) & (
    df.iloc[:, THREE_PRIME_COL].between(115856, 115868)
) & (df.iloc[:, STRAND_COL] == '-/+')

# Condition group 3: four sub-conditions (corrected interval order)
cond3_1 = (
    (df.iloc[:, FIVE_PRIME_COL].between(204252, 204263) |
    df.iloc[:, FIVE_PRIME_COL].between(20365, 204444))  # Corrected to ascending order
) & (
    df.iloc[:, THREE_PRIME_COL].between(207675, 207688)
) & (df.iloc[:, STRAND_COL] == '+')

cond3_2 = (
    df.iloc[:, FIVE_PRIME_COL].between(260980, 260990)
) & (
    df.iloc[:, THREE_PRIME_COL].between(259345, 259364)
) & (df.iloc[:, STRAND_COL] == '-')

cond3_3 = (
    (df.iloc[:, FIVE_PRIME_COL].between(204252, 204263) |
    df.iloc[:, FIVE_PRIME_COL].between(204365, 204444))  # Corrected to ascending order
) & (
    df.iloc[:, THREE_PRIME_COL].between(259345, 259364)
) & (df.iloc[:, STRAND_COL] == '+/-')

cond3_4 = (
    df.iloc[:, FIVE_PRIME_COL].between(260980, 260990)
) & (
    df.iloc[:, THREE_PRIME_COL].between(207675, 207688)
) & (df.iloc[:, STRAND_COL] == '-/+')

# Combine all conditions
all_conditions = (
    cond1_1 | cond1_2 | cond1_3 | cond1_4 |
    cond2_1 | cond2_2 | cond2_3 | cond2_4 |
    cond3_1 | cond3_2 | cond3_3 | cond3_4
)

# Apply filtering conditions
filtered_df = df[all_conditions]

# Save results to new CSV file
filtered_df.to_csv("intron-FL.csv", index=False)

# Output statistics
print(f"Filtering completed! Number of records matching criteria: {len(filtered_df)}")
print(f"Results saved to: intron-FL.csv")