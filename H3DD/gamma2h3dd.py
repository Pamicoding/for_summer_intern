import logging
from pathlib import Path
from datetime import datetime
import csv
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def convert_utc_datetime(utc_datetime_str):
    """
    Converts a UTC datetime string into a datetime object.
    Modify this function based on the format of your UTCDateTime strings.
    """
    return datetime.strptime(utc_datetime_str, '%Y-%m-%dT%H:%M:%S.%f')

def gamma_reorder(ori_csv, reorder_csv):
    # Open the input CSV file and read its contents
    with open(ori_csv, 'r') as file:
        reader = csv.reader(file)
        header = next(reader)
        data = list(reader)
    data.sort(key=lambda row: convert_utc_datetime(row[0]))
    # Write the sorted data to a new CSV file
    with open(reorder_csv, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header) 
        writer.writerows(data) 

def gamma_chunk_split(split_dir, reorder_csv):
    for i, chunk in enumerate(pd.read_csv(reorder_csv, chunksize=4000)): # event needs to < 5000 for executing h3dd
        chunk.to_csv(split_dir / f'gamma_events_{i}.csv', index=False)

def transform(index):
    gamma_picks = picks
    gamma_events = split_dir / f'gamma_events_{index}.csv'
    logging.info(f'we are in gamma_events_{index}') 
    with open(gamma_events,'r') as f:
        lines = f.readlines()
    with open(gamma_picks,'r') as picks_read:
        picks_lines = picks_read.readlines()
    output_file = output_dir / f'gamma_new_{index}.dat_ch'
    with open(output_file,'a') as r:
        for line in lines[1:]:
            item = line.split(',')
            utc_time = datetime.strptime(item[0], '%Y-%m-%dT%H:%M:%S.%f')
            ymd = utc_time.strftime('%Y%m%d')
            hh = utc_time.hour
            mm = utc_time.minute
            ss = round(utc_time.second + utc_time.microsecond / 1000000, 2)
            lon_int = int(float(item[-3]))
            lon_deg = (float(item[-3]) - lon_int)*60
            lat_int = int(float(item[-2]))
            lat_deg = (float(item[-2]) - lat_int)*60
            depth = round(float(item[-1]),2)
            event_index = item[9]
            r.write(f"{ymd:>9}{hh:>2}{mm:>2}{ss:>6.2f}{lat_int:2}{lat_deg:0>5.2f}{lon_int:3}{lon_deg:0>5.2f}{depth:>6.2f}\n")
            for p_line in picks_lines[1:]:
                picks_index = p_line.split(',')[-2]
                if event_index == picks_index:
                    part = p_line.split(',')
                    phase = part[5]
                    sta = part[0].split('.')[1]
                    pick_time = datetime.strptime(part[3], '%Y-%m-%dT%H:%M:%S.%f')
                    if mm == 59 and pick_time.minute == 0: # modify
                        wmm = int(60)
                    else:
                        wmm = pick_time.minute
                    wss = round(pick_time.second + pick_time.microsecond / 1000000, 2)
                    wei = '1.00'
                    if phase == 'P':
                        r.write(f"{' ':1}{sta:<4}{'0.0':>6}{'0':>4}{'0':>4}{wmm:>4}{wss:>6.2f}{'0.01':>5}{wei:>5}{'0.00':>6}{'0.00':>5}{'0.00':>5}\n")
                    else:
                        r.write(f"{' ':1}{sta:<4}{'0.0':>6}{'0':>4}{'0':>4}{wmm:>4}{'0.00':>6}{'0.00':>5}{'0.00':>5}{wss:>6.2f}{'0.01':>5}{wei:>5}\n")
    logging.info(f'gamma_event_{index} transform is done')                    

if __name__ == '__main__':
    gamma_result_dir = Path("/home/patrick/Work/AutoQuake/GaMMA/results/Hualien_0428")
    ori_csv = gamma_result_dir / 'gamma_events.csv'
    reorder_csv = gamma_result_dir / 'gamma_reorder_events.csv'
    split_dir = gamma_result_dir / 'split_dir'
    split_dir.mkdir(parents=True, exist_ok=True)
    picks = gamma_result_dir / 'gamma_picks.csv'
    output_dir = gamma_result_dir / 'for_h3dd'
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(filename='trans.log',level=logging.INFO,filemode='w')

    # reorder
    gamma_reorder(ori_csv, reorder_csv)
    # split
    gamma_chunk_split(split_dir, reorder_csv)
    
    # gamma2h3dd
    chunk_num = sum(1 for _ in split_dir.iterdir() if _.is_file())
    print(f"chunk num: {chunk_num}")
    index_list = np.arange(0, chunk_num)

    with ThreadPoolExecutor(max_workers=chunk_num) as executor:
        executor.map(transform, index_list)
    
