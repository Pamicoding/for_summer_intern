import os
import logging
from pathlib import Path
from datetime import datetime
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
    df = pd.read_csv(ori_csv)
    df['time'] = pd.to_datetime(df['time'])
    df_sort = df.sort_values(by='time')
    df_sort = df_sort.reset_index(drop=True)
    df_sort.to_csv(reorder_csv)

def gamma_chunk_split(split_dir, reorder_csv):
    split_dir.mkdir(parents=True, exist_ok=True)
    for i, chunk in enumerate(pd.read_csv(reorder_csv, chunksize=4000)):
        chunk.to_csv(split_dir / f'gamma_events_{i}.csv',index=False)

def transform(index, split_dir, gamma_picks, output_dir):
    gamma_events = split_dir / f'gamma_events_{index}.csv'
    logging.info(f'we are in gamma_events_{index}') 
    df_catalog = pd.read_csv(gamma_events)
    
    df_catalog['time'] = pd.to_datetime(df_catalog['time'])
    df_catalog['ymd'] = df_catalog['time'].dt.strftime('%Y%m%d')
    df_catalog['hour'] = df_catalog['time'].dt.hour
    df_catalog['minute'] = df_catalog['time'].dt.minute
    df_catalog['seconds'] = df_catalog['time'].dt.second + df_catalog['time'].dt.microsecond / 1_000_000
    df_catalog['lon_int'] = df_catalog['longitude'].apply(lambda x: int(x))
    df_catalog['lon_deg'] = (df_catalog['longitude'].apply(lambda x: float(x)) - df_catalog['lon_int'])*60
    df_catalog['lat_int'] = df_catalog['latitude'].apply(lambda x: int(x))
    df_catalog['lat_deg'] = (df_catalog['latitude'].apply(lambda x: float(x)) - df_catalog['lat_int'])*60
    df_catalog['depth'] = df_catalog['depth_km'].round(2)
    
    df_picks = pd.read_csv(gamma_picks)
    df_picks['phase_time'] = pd.to_datetime(df_picks['phase_time'])
    df_picks['minute'] = df_picks['phase_time'].dt.minute
    df_picks['seconds'] = df_picks['phase_time'].dt.second + df_picks['phase_time'].dt.microsecond / 1_000_000
    
    output_file = output_dir / f'gamma_new_{index}.dat_ch'
    with open(output_file,'w') as r:
        for _, row in df_catalog.iterrows():
            event_index = row['event_index']
            logging.info(f"***Event {event_index}***")
            r.write(f"{row['ymd']:>9}{row['hour']:>2}{row['minute']:>2}{row['seconds']:>6.2f}{row['lat_int']:2}{row['lat_deg']:0>5.2f}{row['lon_int']:3}{row['lon_deg']:0>5.2f}{row['depth']:>6.2f}\n")
            for _, pick_row in df_picks[df_picks['event_index'] == event_index].iterrows():
                if row['minute'] == 59 and pick_row['minute'] == 0: 
                    wmm = int(60)
                else:
                    wmm = pick_row['minute']
                weight = '1.00'
                if pick_row['phase_type'] == 'P':
                    r.write(f"{' ':1}{pick_row['station_id']:<4}{'0.0':>6}{'0':>4}{'0':>4}{wmm:>4}{pick_row['seconds']:>6.2f}{'0.01':>5}{weight:>5}{'0.00':>6}{'0.00':>5}{'0.00':>5}\n")
                else:
                    r.write(f"{' ':1}{pick_row['station_id']:<4}{'0.0':>6}{'0':>4}{'0':>4}{wmm:>4}{'0.00':>6}{'0.00':>5}{'0.00':>5}{pick_row['seconds']:>6.2f}{'0.01':>5}{weight:>5}\n")
                logging.info('next')
    logging.info(f'gamma_event_{index} transform is done')                    

if __name__ == '__main__':
    ori_csv = Path('/home/patrick/Work/EQNet/tests/hualien_0403/result/gamma_events.csv') # gamma_events.csv
    reorder_csv = Path('/home/patrick/Work/EQNet/tests/hualien_0403/gamma_seis_das/gamma_order.csv') # the path of reorder csv you want to name
    split_dir = Path('/home/patrick/Work/EQNet/tests/hualien_0403/gamma_seis_das/split_dir') # the path put the split files
    gamma_picks = Path('/home/patrick/Work/EQNet/tests/hualien_0403/gamma_seis_das/gamma_picks.csv') # gamma_picks.csv
    output_dir = Path('/home/patrick/Work/EQNet/tests/hualien_0403/gamma_seis_das/for_h3dd') # where the output format will generate.

    logging.basicConfig(
        filename='2h3dd.log',
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filemode='w'
        )
    
    gamma_reorder(ori_csv, reorder_csv) # makes gamma event follow the oreder
    gamma_chunk_split(split_dir, reorder_csv) # split the gamma event due to the limitation of h3dd in a execution is 5000 events. (we set 4000) 

    chunk_num = len(os.listdir(split_dir)) # count the num of split files.
    index_list = np.arange(0, chunk_num) # create a list that was utilized to ietrate.
    print("let's go")
    with ThreadPoolExecutor(max_workers=chunk_num) as executor:
        executor.map(transform, index_list, [split_dir] * chunk_num, [gamma_picks] * chunk_num, [output_dir] * chunk_num)

    





