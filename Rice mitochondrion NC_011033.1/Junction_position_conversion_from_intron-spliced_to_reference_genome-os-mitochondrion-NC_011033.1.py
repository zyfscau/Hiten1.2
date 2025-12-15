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
        raw_col = df_original.columns[0]  # Modified: extract column 1 (index 0)
        df_process['1_Original_Column'] = df_original[raw_col]

        
        split_data = df_process['1_Original_Column'].str.split('\|', expand=True)
        split_data.columns = ['MT', 'C', 'D', 'E', 'Strand']

        
        df_process['2_MT'] = split_data['MT']        
        df_process['3_C_str'] = split_data['C']      
        df_process['6_D_str'] = split_data['D']        
        df_process['9_E'] = split_data['E']          
        df_process['10_J'] = split_data['Strand']    

        
        blank_cols = [4,5,7,8,11,12,13,14,15,16,17]
        for pos in blank_cols:
            df_process.insert(pos-1, f'{pos}_Blank', '')

        
        df_process['3_C_num'] = pd.to_numeric(df_process['3_C_str'], errors='coerce').fillna(0)
        df_process['6_D_num'] = pd.to_numeric(df_process['6_D_str'], errors='coerce').fillna(0)

        
        conditions_4 = [
            df_process['3_C_num'] <= 2557,
            (df_process['3_C_num'] > 2557) & (df_process['3_C_num'] <= 7158),
            (df_process['3_C_num'] > 7158) & (df_process['3_C_num'] <= 9932),
            (df_process['3_C_num'] > 9932) & (df_process['3_C_num'] <= 10711),
            (df_process['3_C_num'] > 10711) & (df_process['3_C_num'] <= 17069),
            (df_process['3_C_num'] > 17069) & (df_process['3_C_num'] <= 21824),
            (df_process['3_C_num'] > 21824) & (df_process['3_C_num'] <= 24205),
            (df_process['3_C_num'] > 24205) & (df_process['3_C_num'] <= 28676),
            (df_process['3_C_num'] > 28676) & (df_process['3_C_num'] <= 33976),
        ]

        
        values_4 = [
            9852, 235332, 242480, 305100, 352481, 384717, 290827,
            205314, -10213
        ]

        
        df_process['4_Blank'] = np.select(conditions_4, values_4, default=0).astype(str)

        
        df_process['5_Blank'] = (df_process['3_C_num'] + df_process['4_Blank'].astype(float)).abs()
        df_process['5_Blank'] = df_process['5_Blank'].astype(int).astype(str)

        
        conditions_7 = [
            df_process['6_D_num'] <= 2557,
            (df_process['6_D_num'] > 2557) & (df_process['6_D_num'] <= 7158),
            (df_process['6_D_num'] > 7158) & (df_process['6_D_num'] <= 9932),
            (df_process['6_D_num'] > 9932) & (df_process['6_D_num'] <= 10711),
            (df_process['6_D_num'] > 10711) & (df_process['6_D_num'] <= 17069),
            (df_process['6_D_num'] > 17069) & (df_process['6_D_num'] <= 21824),
            (df_process['6_D_num'] > 21824) & (df_process['6_D_num'] <= 24205),
            (df_process['6_D_num'] > 24205) & (df_process['6_D_num'] <= 28676),
            (df_process['6_D_num'] > 28676) & (df_process['6_D_num'] <= 33976),
        ]

        
        df_process['7_Blank'] = np.select(conditions_7, values_4, default=0).astype(str)
        
        
        df_process['8_Blank'] = (df_process['6_D_num'] + df_process['7_Blank'].astype(float)).abs()
        df_process['8_Blank'] = df_process['8_Blank'].astype(int).astype(str)

        
        def condition_k(x):
            if x <= 2557: return '1'
            elif x <= 10711: return '2'
            elif x <= 17069: return '1'
            elif x <= 28676: return '2'
            elif x <= 33976: return '1'
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
        
        
        df_process['9_E'] = pd.to_numeric(df_process['9_E'], errors='coerce').fillna(0).astype(int).astype(str)
        
        
        df_process['17_Q'] = (
            df_process['2_MT'] + '|' +
            df_process['5_Blank'] + '|' +
            df_process['8_Blank'] + '|' +
            df_process['9_E'] + '|' +
            df_process['15_O']
        )

        
        final_cols = [
            '1_Original_Column', '2_MT', '3_C_str', '4_Blank', '5_Blank', 
            '6_D_str', '7_Blank', '8_Blank', '9_E', '10_J',
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