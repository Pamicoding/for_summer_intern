def append(ori_list, process_list):
    with open(ori_list, 'r') as infile:
        with open(process_list, 'w') as outfile:
            for line in infile:
                new_line = f"{line.strip()} 19010101 21001231\n"
                outfile.write(new_line)


if __name__ == '__main__':
    ori_list = '/home/patrick/Work/EQNet/tests/hualien_0403/station.txt'
    process_list = '/home/patrick/Work/for_summer_intern/H3DD/test.txt'
    append(ori_list, process_list)