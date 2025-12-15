import os
import glob
import pandas as pd
import numpy as np

def process_position_correction():
    
    for input_file in glob.glob('*vlookup-result-sift*'):
        
        output_file = input_file
        
        process_file = f"{os.path.splitext(input_file)[0]}_process.csv"
        
        
        df_original = pd.read_csv(input_file, dtype=str)

        
        df_process = pd.DataFrame(index=df_original.index)
        raw_col = df_original.columns[0]  
        df_process['1_original_column'] = df_original[raw_col]

        
        split_data = df_process['1_original_column'].str.split('\|', expand=True)
        split_data.columns = ['MT', 'C', 'D', 'E', 'Strand']

        df_process['2_MT'] = split_data['MT']
        df_process['3_C_str'] = split_data['C']
        df_process['6_D_str'] = split_data['D']
        df_process['9_E'] = split_data['E']
        df_process['10_J'] = split_data['Strand']

        
        blank_cols = [4,5,7,8,11,12,13,14,15,16,17]
        for pos in blank_cols:
            df_process.insert(pos-1, f'{pos}_blank', '')

        
        df_process['3_C_num'] = pd.to_numeric(df_process['3_C_str'], errors='coerce').fillna(0)
        df_process['6_D_num'] = pd.to_numeric(df_process['6_D_str'], errors='coerce').fillna(0)

        # position correction conditions
        conditions_4 = [
            df_process['3_C_num'] <= 2385,
            (df_process['3_C_num'] > 2385) & (df_process['3_C_num'] <= 2467),
            (df_process['3_C_num'] > 2467) & (df_process['3_C_num'] <= 2658),
            (df_process['3_C_num'] > 2658) & (df_process['3_C_num'] <= 2716),
            (df_process['3_C_num'] > 2716) & (df_process['3_C_num'] <= 3975),
            (df_process['3_C_num'] > 3975) & (df_process['3_C_num'] <= 6128),
            (df_process['3_C_num'] > 6128) & (df_process['3_C_num'] <= 6520),
            (df_process['3_C_num'] > 6520) & (df_process['3_C_num'] <= 6681),
            (df_process['3_C_num'] > 6681) & (df_process['3_C_num'] <= 7254),
            (df_process['3_C_num'] > 7254) & (df_process['3_C_num'] <= 8442),
            (df_process['3_C_num'] > 8442) & (df_process['3_C_num'] <= 10672),
            (df_process['3_C_num'] > 10672) & (df_process['3_C_num'] <= 11888),
            (df_process['3_C_num'] > 11888) & (df_process['3_C_num'] <= 11910),
            (df_process['3_C_num'] > 11910) & (df_process['3_C_num'] <= 12305),
            (df_process['3_C_num'] > 12305) & (df_process['3_C_num'] <= 13455),
        ]

        # position correction value
        values_4 = [
            48489, 318543, 319936, 546056, -269949, -517360, -516538,
            105093, 106998, 108384, 194189, 195059, -278504, -272853, -271902
        ]

        df_process['4_blank'] = np.select(conditions_4, values_4, default=0).astype(str)

        
        df_process['5_blank'] = (df_process['3_C_num'] + df_process['4_blank'].astype(float)).abs()
        df_process['5_blank'] = df_process['5_blank'].astype(int).astype(str)

        
        conditions_7 = [
            df_process['6_D_num'] <= 2385,
            (df_process['6_D_num'] > 2385) & (df_process['6_D_num'] <= 2467),
            (df_process['6_D_num'] > 2467) & (df_process['6_D_num'] <= 2658),
            (df_process['6_D_num'] > 2658) & (df_process['6_D_num'] <= 2716),
            (df_process['6_D_num'] > 2716) & (df_process['6_D_num'] <= 3975),
            (df_process['6_D_num'] > 3975) & (df_process['6_D_num'] <= 6128),
            (df_process['6_D_num'] > 6128) & (df_process['6_D_num'] <= 6520),
            (df_process['6_D_num'] > 6520) & (df_process['6_D_num'] <= 6681),
            (df_process['6_D_num'] > 6681) & (df_process['6_D_num'] <= 7254),
            (df_process['6_D_num'] > 7254) & (df_process['6_D_num'] <= 8442),
            (df_process['6_D_num'] > 8442) & (df_process['6_D_num'] <= 10672),
            (df_process['6_D_num'] > 10672) & (df_process['6_D_num'] <= 11888),
            (df_process['6_D_num'] > 11888) & (df_process['6_D_num'] <= 11910),
            (df_process['6_D_num'] > 11910) & (df_process['6_D_num'] <= 12305),
            (df_process['6_D_num'] > 12305) & (df_process['6_D_num'] <= 13455),
        ]

        df_process['7_blank'] = np.select(conditions_7, values_4, default=0).astype(str)
        
        
        df_process['8_blank'] = (df_process['6_D_num'] + df_process['7_blank'].astype(float)).abs()
        df_process['8_blank'] = df_process['8_blank'].astype(int).astype(str)

        
        def condition_k(x):
            if x <= 2716: return '1'
            elif x <= 6520: return '2'
            elif x <= 11888: return '1'
            elif x <= 13455: return '2'
            else: return '0'
        df_process['11_K'] = df_process['3_C_num'].apply(condition_k)
        df_process['12_L'] = df_process['6_D_num'].apply(condition_k)

        
        df_process = df_process.fillna('')

        
        df_process['10_J'] = df_process['10_J'].astype(str).str.strip().fillna('')
        symbol_map = {
            '+': ['+', '正', 'positive'],
            '-': ['-', '负', 'negative']
        }
        symbol_dict = {v: k for k, values in symbol_map.items() for v in values}
        df_process['10_J'] = df_process['10_J'].map(lambda x: symbol_dict.get(x, x))

        
        df_process['13_M'] = df_process['10_J'] + df_process['11_K'] + df_process['12_L']
        df_process['14_N'] = df_process['13_M']

        
        def final_judge(row):
            j, k, l = row['10_J'], row['11_K'], row['12_L']
            
            if j == '+' and k == '1' and l == '1':
                return '+'
            elif j == '+' and k == '1' and l == '2':
                return '+/-'
            elif j == '+' and k == '2' and l == '2':
                return '-'
            elif j == '+' and k == '2' and l == '1':
                return '-/+'
            elif j == '-' and k == '1' and l == '1':
                return '-'
            elif j == '-' and k == '1' and l == '2':
                return '-/+'
            elif j == '-' and k == '2' and l == '2':
                return '+'
            elif j == '-' and k == '2' and l == '1':
                return '+/-'
            else:
                return ''
        
        df_process['15_O'] = df_process.apply(final_judge, axis=1)

        df_process['16_P'] = '|'
        
        
        df_process['9_E'] = pd.to_numeric(df_process['9_E'], errors='coerce').fillna(0).astype(int).astype(str)
        
        
        df_process['17_Q'] = (
            df_process['2_MT'] + '|' +
            df_process['5_blank'] + '|' +
            df_process['8_blank'] + '|' +
            df_process['9_E'] + '|' +
            df_process['15_O']
        )

        
        final_cols = [
            '1_original_column', '2_MT', '3_C_str', '4_blank', '5_blank', 
            '6_D_str', '7_blank', '8_blank', '9_E', '10_J',
            '11_K', '12_L', '13_M', '14_N', '15_O',
            '16_P', '17_Q'
        ]
        df_process = df_process[final_cols]

        
        df_process.to_csv(process_file, index=False, encoding='utf-8-sig')

        
        df_original.iloc[:, 0] = df_process['17_Q']
        df_original.to_csv(output_file, index=False, encoding='utf-8-sig')
        
        
        if os.path.exists(process_file):
            os.remove(process_file)
            print(f"Deleted temporary file: {process_file}")        
        
        
        print(f"Processed and updated file: {output_file}")           

if __name__ == '__main__':
    process_position_correction()