%% Clean up the tags and logs associated with a job

function [] = sub_clean_job(path_logs,name_job)
files{1} = [path_logs filesep name_job '.log'];
files{2} = [path_logs filesep name_job '.finished'];
files{3} = [path_logs filesep name_job '.failed'];
files{4} = [path_logs filesep name_job '.running'];
files{5} = [path_logs filesep name_job '.exit'];
files{6} = [path_logs filesep name_job '.eqsub'];
files{7} = [path_logs filesep name_job '.oqsub'];
files{8} = [path_logs filesep name_job '.profile.mat'];
files{9} = [path_logs filesep 'tmp' filesep name_job '.sh'];
for num_f = 1:length(files)
    if psom_exist(files{num_f})
        delete(files{num_f});
    end
end
