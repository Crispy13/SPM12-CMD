import os,re,glob,sys,argparse,tempfile
from subprocess import Popen, PIPE, run
import datetime,shlex
import nibabel, nibabel.processing
import gzip, shutil

import sys
if sys.version_info[0] < 3:
    raise Exception("Python 3.0+ is needed.")

# Get Arguments
parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="A script to execute spm12 at command line.")

# Positional arugment : command
subparsers=parser.add_subparsers(metavar= 'command', required = True)
   
parser_convert = subparsers.add_parser('convert', help='convert .m script into commandline-available script.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_convert.set_defaults(command = 'convert')
parser_convert.add_argument('m_script', type=str, help='m script file to convert.')
parser_convert.add_argument('output', type=str, help='output file path.')

parser_gunzip = subparsers.add_parser('gunzip', help='gunzip nii.gz files.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_gunzip.set_defaults(command = 'gunzip')
parser_gunzip.add_argument("root_folder", help="root folder to search mri files.")
parser_gunzip.add_argument("patient_id", default = 'all', help="Patient ID.", nargs='*')
parser_gunzip.add_argument("-s", "--scan-types", default = ['t1ce', 't2', 'seg'], 
                        help="A list of available image scan types.(e.g. -s t1ce t2 seg) The last element of this should be 'seg'. Otherwise it will malfunction.",
                        type=str, nargs='+')

parser_resample = subparsers.add_parser('resample', help='resample a mri image.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_resample.set_defaults(command = 'resample')
parser_resample.add_argument('input_path', type=str, help='nii file to resample.')
parser_resample.add_argument('-vs', "--voxel_size", default=[1,1,1], type=int, help='voxel size(e.g. -vs 2 2 2)', nargs='+')

parser_run = subparsers.add_parser('run', help='run spm12 process pipeline.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_run.set_defaults(command = 'run')
parser_run.add_argument("root_folder", help="root folder to search mri files.")
parser_run.add_argument("patient_id", default = 'all', help="Patient ID.", nargs='*')
parser_run.add_argument("-s", "--scan-types", default = ['t1ce', 't2', 'seg'], 
                        help="A list of available image scan types.(e.g. -s t1ce t2 seg) The last element of this should be 'seg'. Otherwise it will malfunction.",
                        type=str, nargs='+')
parser_run.add_argument("-m", "--m-script-archive", default = "/data/eck/CustomScriptsArchive/spm12_scripts", 
                        help="m-script-archive.",
                        type=str, nargs='?')

parser.add_argument("--mcr-path", default = "/data/eck/software/MCR/v97", 
                        help="m-script-archive.",
                        type=str, nargs='?')
parser.add_argument("-l", "--log-file", default = "spm12_cmdline.log", 
                    help="log file.",
                    type=str, nargs='?')

# parser_b = subparsers.add_parser('b', help='b help')
# parser_b.add_argument('--baz', choices='XYZ', help='baz help')

args=parser.parse_args()

### Create Logger

from crispy13 import create_logger

logger = create_logger(file_handler_path = args.log_file)

from crispy13.core.ecf import cus_excepthook

sys.excepthook=cus_excepthook(logger)

### Input dict functions

def if_found_return_groups(pattern, iterable, group_index = None, flags = 0):

    r = []
    for i in iterable:
        sr1 = re.search(pattern, i, flags = flags) #search result 
        
        if sr1:
            r.append(sr1.groups()) if group_index is None else r.append(sr1.group(group_index))
        else:
            pass
        
    return r

class spm12_input_dict:
    def __init__(self, root_folder, img_types, pid = None, file_extension = ".nii"):
        self.root_folder = os.path.abspath(root_folder)
        self.img_types = img_types
        self.file_extension = file_extension
        self.fep = re.sub("\.", "\\.", file_extension) # file_extension_pattern
        self.set_input_dict(root_folder, img_types, pid)
        
        
    def __str__(self):
        return self.values.__str__()
    
    def __repr__(self):
        return self.values.__repr__()
    
    def set_input_dict(self, root_folder, img_types, pid):
        root_folder = os.path.abspath(root_folder)
        tfl1 = glob.glob(f"{root_folder}/**/*{self.file_extension}", recursive = True) # temp folder list 1
        
        # modify pid variable.
        if (pid == 'all') or (pid == ['all']):
            pid = set(if_found_return_groups("[0-9]{8}", tfl1, 0))
        elif isinstance(pid, int):
            pid = [str(pid)]
        elif pid is None:
            raise ValueError(f"pid should be given.")

        self.pid = pid
        input_dict = dict.fromkeys(self.pid)
        
        for i in input_dict:
            gr1 = glob.glob(f"{root_folder}/**/{i}/**/*{self.file_extension}", recursive = True)
            if len(gr1) == 0: raise ValueError(f"gr1 has nothing. {i}\n{gr1}")
            
            if self.file_extension != "_final.nii":
                gr1 = list(filter(lambda x:re.search("/(?:[rmc]{1,2}|mean).*" + self.fep, x) is None, gr1))
                gr1 = list(filter(lambda x:re.search("/.*_final" + self.fep, x) is None, gr1))
            else:
                pass
            
            if len(gr1) == 0: raise ValueError(f"gr1_2 has nothing. {i}\n{gr1}")
            
            ifr1 = set(if_found_return_groups("([0-9]{4}-[0-9]{2}-[0-9]{2})", gr1, group_index = 1))
            if len(ifr1) == 0: raise ValueError(f"ifr1 has nothing. {i}\nifr1:{ifr1}\ngr1:{gr1}")
            
            input_dict[i] = dict.fromkeys(ifr1)
            
            for j in input_dict[i]:
                td = {}
                    
                for t in img_types:
                    if t != 'seg':
                        tt= 't1' if t=='t1ce' else t
                        tl = if_found_return_groups("^((?!roi|seg|label).)*$",
                                                        if_found_return_groups(f"{root_folder}(?:/[^/\n]+)*/{i}/{j}.*/[^/\n]*{tt}[^/\n]*{self.fep}", gr1, 0, re.I),
                                                        0,
                                                        re.I)
                        if len(tl) == 1:
                            td[t] = tl[0]
                        else:
                            raise ValueError(f"More than one path was found. {tl} {i} {j} {t} {tt} {gr1}")

                    else:
                        tl = if_found_return_groups(f"{root_folder}(?:/[^/\n]+)*/{i}/{j}.*/[^/\n]*(?:roi|seg|label)[^/\n]*{self.fep}", gr1, 0, re.I)
                        if len(tl) == 1:
                            td[t] = tl[0]
                        else:
                            raise ValueError(f"More than one path was found. {tl} {i} {j} {t} {tt} {gr1}")

                input_dict[i][j] = td
            
        self.values = input_dict

### Gunzip functions

def flatten_dict(d):
    """
    https://stackoverflow.com/questions/52081545/python-3-flattening-nested-dictionaries-and-lists-within-dictionaries
    """
    
    out = {}
    for key, val in d.items():
        if isinstance(val, dict):
            val = [val]
        if isinstance(val, list):
            for subdict in val:
                deeper = flatten_dict(subdict).items()
                out.update({key + '_' + key2: val2 for key2, val2 in deeper})
        else:
            out[key] = val
    return out

def gunzip_nii_gz_files(root_folder, scan_types, patient_id):
    nii_gz_dict = spm12_input_dict(root_folder, scan_types, patient_id, ".nii.gz")
    nii_gz_dict.values = flatten_dict(nii_gz_dict.values)
    
    for f in nii_gz_dict.values.values():
        nf = re.sub("\.nii\.gz", ".nii", f)
        logger.info(f"Gunzip {f} ...")
        with gzip.open(f, 'r') as f_in, open(nf, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
    logger.info(f">>> Gunzip complete.")

def gzip_nii_files(root_folder, scan_types, patient_id, nii_file_pattern = "_final.nii"):
    nii_dict = spm12_input_dict(root_folder, scan_types, patient_id, nii_file_pattern)
    nii_dict.values = flatten_dict(nii_dict.values)
    
    for f in nii_dict.values.values():
        nf = re.sub("\.nii", ".nii.gz", f)
        logger.info(f"Gzip {f} ...")
        with open(f, 'rb') as f_in, gzip.open(nf, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
    logger.info(f">>> Gzip complete.")

### Convert Functions

def convert_mscript(process, inputs):
    inputs = list(map(lambda x:f"'{x}'", inputs))
    with open(f"{args.m_script_archive}/{process}.m", 'r') as f:
            contents = f.read()
            
    logger.info(f"Before converting script:\n{contents}\n")
    
    # Edit Realignment script
    # inputs will be T1 + seg or T2
    if process == 'realignment':
        inputs_string = '\n' + '\n'.join(inputs) + '\n'
        ir = re.sub("\{[^{]*'.+'[^}]*\}", f"{{{inputs_string}}}", contents)
    
    # Edit Coregistration script
    # inputs will be T1 + T2.
    if process == 'coregistration':
        ir = re.sub("(\.ref = ){'.+'}", f"\g<1>{{{inputs[0]}}}", contents)
        ir = re.sub("(\.source = ){'.+'}", f"\g<1>{{{inputs[1]}}}", ir)
        try:
            input_strings = '\n' + '\n'.join(inputs[2:]) + '\n'
            ir = re.sub("(\.other = ){'.*'}", f"\g<1>{{{input_strings}}}", ir)
        except IndexError:
            pass
        
        
    # Edit segment script
    # inputs will be T1 + T2.
    if process == 'segment':
        ir = re.sub("(\.channel\(1\).vols = ){'.+'}", f"\g<1>{{{inputs[0]}}}", contents)
        ir = re.sub("(\.channel\(2\).vols = ){'.+'}", f"\g<1>{{{inputs[1]}}}", ir)
        
        
    # Edit imgcalc script
    # inputs will be T1 + T2.
    if process == 'imgcalc':
        inputs_string = '\n' + '\n'.join(inputs[:4]) + '\n'

        ir = re.sub("(\.imcalc.input = )\{(?:[^{]*'.+'[^{]*)*\}", f"\g<1>{{{inputs_string}}}", contents)
        ir = re.sub("(\.imcalc.output = )'.+'", f"\g<1>{inputs[4]}", ir)
        ir = re.sub("(\.imcalc.outdir = ){'.+'}", f"\g<1>{{{inputs[5]}}}", ir)
        
    edited_contents = f'''
spm defaults fmri
spm_jobman initcfg

{ir}
spm_jobman('run',matlabbatch);
    '''
    
    logger.info(f"After converting script:\n{edited_contents}\n")
    logger.info(f'''Running {process.upper()} with inputs:{inputs}''')
    
    ts_file = tempfile.NamedTemporaryFile(prefix = f"{process}_") # temp script file
    tfn = ts_file.name
    ts_file.close()
    
    with open(f"{tfn}.m", 'w') as f:
        f.write(edited_contents)
        
    return f"{tfn}.m"

### Resample functions

def resample_mri(input_path, voxel_size = [1, 1, 1], output_path_modifier = "/\g<1>_resampled\g<2>"):
    """
    Resample mri nii file and save it to a file.
    """
    
    output_path = re.sub("/(.+)(\.[^.]+?)$", output_path_modifier, input_path)
    logger.info(f"output_path: {output_path}")
    
    input_img = nibabel.load(input_path)
    resampled_img = nibabel.processing.resample_to_output(input_img, voxel_size)
    nibabel.save(resampled_img, output_path)

### Run functions

def _run_m_script(m_script, MCR_path=args.mcr_path):
    cmd1 = f'''run_spm12.sh {MCR_path} script {m_script}'''
    logger.info(f"running the following command:\n{cmd1}")
          
    p1 = Popen(shlex.split(cmd1), stdout=PIPE)
    p1_stdout, p1_stderr = p1.communicate()
    
    if p1.returncode != 0:
        logger.info(f"Stdout:\n{p1_stdout}")
        logger.info(f"Stderr:\n{p1_stderr}")
        logger.info(f"Return code: {p1.returncode}")
        raise Exception(f"Job failed.")

def run_spm12_process(input_dict):
    if input_dict.__class__.__name__ != spm12_input_dict.__name__: raise ValueError("input_dict should be an instance of spm12_input_dict class.")
    
    pid = input_dict.pid
    input_data = input_dict.values
    img_types = input_dict.img_types
    
    for p in pid:
        logger.info(f"Current patient id: {p}.")
        
        for ts in input_data[p]: # ts = time series
            # run realignment
#             for t in img_types[:-1]:
#                 if t == 't1ce':
#                     inputs = [input_data[p][ts][t]]
#                     m_script = convert_mscript("realignment", inputs)
#                     _run_m_script(m_script) # T1 & label realignment
#                     os.remove(m_script)
#                 else:
#                     inputs = [input_data[p][ts][t]]
#                     m_script = convert_mscript("realignment", inputs)
#                     _run_m_script(m_script) # T2 realignment
#                     os.remove(m_script)
                    
            # run coregistration
            inputs = [input_data[p][ts][t] for t in img_types[:]]
            m_script = convert_mscript("coregistration", inputs)
            _run_m_script(m_script) # T1 & T2 coregistration
            os.remove(m_script)
            
            shutil.move(inputs[-1], re.sub("/([^\n/]*)\.nii", "/\g<1>_final.nii", inputs[-1])) # For consistency of the names of final products.
        
            # run segment
            # Make inputs
            inputs = []
            for i in range(len(img_types[:-1])):
                t = img_types[i]
                if i < 1:
                    inputs.append(input_data[p][ts][t])
                else:
                    inputs.append(re.sub(f"/([^\n/]+\.nii)", "/r\g<1>", input_data[p][ts][t]))
                
            m_script = convert_mscript("segment", inputs)
            _run_m_script(m_script) # T1 & T2 segment
            os.remove(m_script)
            
            # run imgcalc(skull stripping)
            for i in range(len(img_types[:-1])):
                t = img_types[i]
                if i < 1:
                    inputs = re.sub(f"/([^\n/]+\.nii)", "/m\g<1>", input_data[p][ts][t])
                else:
                    inputs = re.sub(f"/([^\n/]+\.nii)", "/mr\g<1>", input_data[p][ts][t])
                    
                inputs = [inputs] + [re.sub(f"/([^\n/]+\.nii)", "/c1\g<1>", input_data[p][ts]['t1ce'])]\
                                + [re.sub(f"/([^\n/]+\.nii)", "/c2\g<1>", input_data[p][ts]['t1ce'])]\
                                + [re.sub(f"/([^\n/]+\.nii)", "/c3\g<1>", input_data[p][ts]['t1ce'])]\
                                + [re.sub(f"^.*/([^\n/]+)\.nii", "\g<1>_final.nii", inputs)]\
                                + [re.search(f"(^.+)/([^\n/]+\.nii)", inputs).group(1)]
            
                m_script = convert_mscript("imgcalc", inputs)
                _run_m_script(m_script) # T1 & T2 imgcalc
                os.remove(m_script)
            
    logger.info(f">> SPM Process complete.")

def run_spm12(root_folder, img_types, pid):
    img_type_ref = ['t1ce', 't1', 't2', 'flair', 'seg']
    for i in img_types:
        if i not in img_type_ref: raise ValueError(f"Invaild img type. Ref = {img_type_ref}.")
    
    if 'seg' in img_types:
        img_types.remove('seg')
        img_types.insert(len(img_types), 'seg')
    if 't1ce' in img_types:
        img_types.remove('t1ce')
        img_types.insert(0, 't1ce')
    
    gunzip_nii_gz_files(root_folder, img_types, pid)
    input_dict = spm12_input_dict(root_folder, img_types, pid)
    run_spm12_process(input_dict)
    gzip_nii_files(root_folder, img_types, pid, nii_file_pattern = "_final.nii")
    
    logger.info(f">> The whole processes complete.")

if args.command == 'convert':
    convert_mscript(args.m_script, args.output)
    logger.info("Converting process complete.")

elif args.command == 'gunzip':
    gunzip_nii_gz_files(args.root_folder, args.scan_types, args.patient_id)

elif args.command == 'resample':
    resample_mri(args.input_path, args.voxel_size)

elif args.command == 'run':
    logger.info(f"patient id : {args.patient_id}")
    subpcs = []
    args.root_folder = os.path.abspath(args.root_folder)
    args.m_script_archive = os.path.abspath(args.m_script_archive)
    
    
    if (len(args.patient_id) > 1) or (args.patient_id == ['all']):
        pids = set(if_found_return_groups("[0-9]{8}", glob.glob(f"{args.root_folder}/**/*.nii.gz", recursive = True), 0))
        for pid in pids:
            print(f'''python {args.m_script_archive}/spm12_cmdline_main.py -l {args.root_folder}/logs/{pid}_spm12_run.log run {args.root_folder} {pid} -s {' '.join(args.scan_types)}'''.split())
            subpcs.append(
                Popen(f'''python {args.m_script_archive}/spm12_cmdline_main.py -l {args.root_folder}/logs/{pid}_spm12_run.log run {args.root_folder} {pid} -s {' '.join(args.scan_types)}'''.split(), stdout = PIPE, stderr = PIPE, cwd = os.getcwd())
            )

        exit_codes =[]
        for sp in subpcs:
            sp.communicate()
            exit_codes.append(sp.returncode)
            
        logger.info(f"Job orders: {pids}")
        logger.info(f"Exit codes: {exit_codes}")
    
    elif (len(args.patient_id) == 1) and (re.search("[0-9]{8}", args.patient_id[0]) is not None):
        run_spm12(args.root_folder, args.scan_types, args.patient_id)

