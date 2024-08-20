# make sure all required imports are done in Snakefile

def convertToMb(string):
    '''
    This function can convert text in the form
    xxG to mb
    If it does not end with G, it returns the string
    It does not handle other cases of invalid input
    '''
    if string.endswith('G'):
        number = int(string.split('G')[0])
        return(number*1000)
    else:
        return(string)

def convertToSec(string):
    '''
    This function can convert text in the form
    D-hh:mm:ss to seconds
    D - # days
    hh # hours
    mm # mins
    ss # secs
    '''
    days = string.split('-')[0]
    hrs = string.split('-')[1].split(':')[0]
    min = string.split('-')[1].split(':')[1]
    sec = string.split('-')[1].split(':')[2]
    total = int(sec)
    total = total + 60*int(min)
    total = total + 60*60*int(hrs)
    total = total + 24*60*60*int(days)
    return(total)

def remove_run_name_fluff(string):
    '''
    This function removes the run name fluff from the file name
    '''
    start = string.find('EKDN')
    end = string.find('SXC_')
    return(string[:start] + string[end+4:])

raw_paths_dict = {x: glob.glob(f'results/00_RawData/from_Novogene/{x}/*.fq.gz') for x in SAMPLES}
name_clean_dict = {remove_run_name_fluff(x): x for x in list(chain(*raw_paths_dict.values()))}

def get_input_file(sample, lane, read):
    '''
    This function returns the path of the file given the sample, read, lane and run
    '''
    the_path = name_clean_dict[f'results/00_RawData/from_Novogene/{sample}/{sample}_{lane}_{read}.fq.gz']
    return(the_path)

trimmed_paths_dict = {x: [os.path.join(f'results/01_TrimmingFiltering', remove_run_name_fluff(os.path.basename(y))) for y in glob.glob(f'results/00_RawData/from_Novogene/{x}/*.fq.gz')] for x in SAMPLES}

def get_trimmed_files(sample, read):
    '''
    This function returns the path of all the trimmed files given the sample and read
    since there are multiple lanes, it returns a list that can be used in the rule
    running bowtie which asks for them to be a comma separated list
    '''
    the_files = [x for x in trimmed_paths_dict[sample] if x.endswith(f'_{read}.fq.gz')]
    return(the_files)

    