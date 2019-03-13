def split_xl(input_list):
    split_point=0    
    for i in range(len(input_list)):
        if i<len(input_list)-1:
            if input_list[i]=='-' and input_list[i+1]=='-':
                split_point = i;
    return split_point