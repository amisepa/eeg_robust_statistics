function [] = sub_add_line_log(file_write,str_write,flag_verbose)
if flag_verbose
    fprintf('%s',str_write)
end
if ischar(file_write)
    hf = fopen(file_write,'a');
    fprintf(hf,'%s',str_write);
    fclose(hf);
else
    fprintf(file_write,'%s',str_write);
end
