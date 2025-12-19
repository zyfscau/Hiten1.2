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

        #  position correction conditions
        conditions_4 = [
            df_process['3_C_num'] <= 2114,
            (df_process['3_C_num'] > 2114) & (df_process['3_C_num'] <= 2346),
            (df_process['3_C_num'] > 2346) & (df_process['3_C_num'] <= 3203),
            (df_process['3_C_num'] > 3203) & (df_process['3_C_num'] <= 3980),
            (df_process['3_C_num'] > 3980) & (df_process['3_C_num'] <= 5736),
            (df_process['3_C_num'] > 5736) & (df_process['3_C_num'] <= 7776),
            (df_process['3_C_num'] > 7776) & (df_process['3_C_num'] <= 8994),
            (df_process['3_C_num'] > 8994) & (df_process['3_C_num'] <= 11139),
            (df_process['3_C_num'] > 11139) & (df_process['3_C_num'] <= 12546),
            (df_process['3_C_num'] > 12546) & (df_process['3_C_num'] <= 14670),
            (df_process['3_C_num'] > 14670) & (df_process['3_C_num'] <= 14900),
            (df_process['3_C_num'] > 14900) & (df_process['3_C_num'] <= 16059),
            (df_process['3_C_num'] > 16059) & (df_process['3_C_num'] <= 18065),
            (df_process['3_C_num'] > 18065) & (df_process['3_C_num'] <= 18955),
            (df_process['3_C_num'] > 18955) & (df_process['3_C_num'] <= 20429),
            (df_process['3_C_num'] > 20429) & (df_process['3_C_num'] <= 22438),
            (df_process['3_C_num'] > 22438) & (df_process['3_C_num'] <= 23840),
            (df_process['3_C_num'] > 23840) & (df_process['3_C_num'] <= 26230),
            (df_process['3_C_num'] > 26230) & (df_process['3_C_num'] <= 27662),
            (df_process['3_C_num'] > 27662) & (df_process['3_C_num'] <= 30211),
            (df_process['3_C_num'] > 30211) & (df_process['3_C_num'] <= 31751)
        ]

        # position correction value
        values_4 = [
            -71428, -95241, -94701, -94701, -93397, -13332, -12476,
            24093, 24923, -60725, -59989, -59258, 55087, 55787, 56531,
            -103421, -102368, -110278, -109615, -146162, 145139
        ]

        
        df_process['4_blank'] = np.select(conditions_4, values_4, default=0).astype(str)

        # Modification point 1: Convert column 5_blank to integer string
        df_process['5_blank'] = (df_process['3_C_num'] + df_process['4_blank'].astype(float)).abs()
        df_process['5_blank'] = df_process['5_blank'].astype(int).astype(str)

        
        conditions_7 = [
            df_process['6_D_num'] <= 2114,
            (df_process['6_D_num'] > 2114) & (df_process['6_D_num'] <= 2346),
            (df_process['6_D_num'] > 2346) & (df_process['6_D_num'] <= 3203),
            (df_process['6_D_num'] > 3203) & (df_process['6_D_num'] <= 3980),
            (df_process['6_D_num'] > 3980) & (df_process['6_D_num'] <= 5736),
            (df_process['6_D_num'] > 5736) & (df_process['6_D_num'] <= 7776),
            (df_process['6_D_num'] > 7776) & (df_process['6_D_num'] <= 8994),
            (df_process['6_D_num'] > 8994) & (df_process['6_D_num'] <= 11139),
            (df_process['6_D_num'] > 11139) & (df_process['6_D_num'] <= 12546),
            (df_process['6_D_num'] > 12546) & (df_process['6_D_num'] <= 14670),
            (df_process['6_D_num'] > 14670) & (df_process['6_D_num'] <= 14900),
            (df_process['6_D_num'] > 14900) & (df_process['6_D_num'] <= 16059),
            (df_process['6_D_num'] > 16059) & (df_process['6_D_num'] <= 18065),
            (df_process['6_D_num'] > 18065) & (df_process['6_D_num'] <= 18955),
            (df_process['6_D_num'] > 18955) & (df_process['6_D_num'] <= 20429),
            (df_process['6_D_num'] > 20429) & (df_process['6_D_num'] <= 22438),
            (df_process['6_D_num'] > 22438) & (df_process['6_D_num'] <= 23840),
            (df_process['6_D_num'] > 23840) & (df_process['6_D_num'] <= 26230),
            (df_process['6_D_num'] > 26230) & (df_process['6_D_num'] <= 27662),
            (df_process['6_D_num'] > 27662) & (df_process['6_D_num'] <= 30211),
            (df_process['6_D_num'] > 30211) & (df_process['6_D_num'] <= 31751)
        ]

        
        df_process['7_blank'] = np.select(conditions_7, values_4, default=0).astype(str)
        
        # Modification point 2: Convert column 8_blank to integer string
        df_process['8_blank'] = (df_process['6_D_num'] + df_process['7_blank'].astype(float)).abs()
        df_process['8_blank'] = df_process['8_blank'].astype(int).astype(str)  # Added integer conversion

        
        def condition_k(x):
            if x <= 8994: return '2'
            elif x <= 12546: return '1'
            elif x <= 16059: return '2'
            elif x <= 20429: return '1'
            elif x <= 31751: return '2'
            else: return '0'
        df_process['11_K'] = df_process['3_C_num'].apply(condition_k)
        df_process['12_L'] = df_process['6_D_num'].apply(condition_k)

        
        df_process = df_process.fillna('')

        
        df_process['10_J'] = df_process['10_J'].astype(str).str.strip().fillna('')
        symbol_map = {
            '+': ['+', 'positive'],
            '-': ['-', 'negative']
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
        
        # Modification point 3: Process numeric format for column 9_E
        df_process['9_E'] = pd.to_numeric(df_process['9_E'], errors='coerce').fillna(0).astype(int).astype(str)
        
        # Modification point 4: Final concatenation
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