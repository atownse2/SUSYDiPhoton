import os

import re

condor_dir = f"/scratch365/atownse2/SUSYDiPhoton/condor"

def submit_batch(executable, arguments, job_name):
    import htcondor
    skim_events = htcondor.Submit({
        "executable": executable,
        "arguments": arguments,
        "output": f"{condor_dir}/out/{job_name}.out",
        "error" : f"{condor_dir}/err/{job_name}.err",
        "log"   : f"{condor_dir}/log/{job_name}.log",              
        "request_cpus": "1",
        "request_memory": "512MB",
        "request_disk": "256MB",
    })
    print(f"Submitting job {job_name}")
    schedd = htcondor.Schedd()
    submit_result = schedd.submit(skim_events)

def submit_local(executable, arguments, job_name):
    command = f'{executable} {arguments}'
    print(command)
    os.system(command)

def clear_logs(job_name):
    for tag in ["out", "err", "log"]:
        if os.path.exists(f"{condor_dir}/{tag}/{job_name}.{tag}"):
            os.remove(f"{condor_dir}/{tag}/{job_name}.{tag}")


def check_logs(job_tag):
    log_files = [f for f in os.listdir(f"{condor_dir}/log") if job_tag in f]
    # err_files = [f for f in os.listdir(f"{condor_dir}/err") if job_tag in f]

    no_exit_code = 0
    n_bad_jobs = 0
    n_total_jobs = len(log_files)
    exit_codes = {}

    for log_file in log_files:
        with open(f"{condor_dir}/log/{log_file}") as f:
            txt = f.read()
            match = re.search(r'with exit-code (\d+)', txt)
            if not match:
                no_exit_code += 1
                continue
            exit_code = int(match.group(1))
            if exit_code != 0:
                if exit_code not in exit_codes:
                    exit_codes[exit_code] = {
                        "file" :f"{condor_dir}/err/{log_file.replace('log', 'err')}",
                        "count" : 1}
                else:
                    exit_codes[exit_code]["count"] += 1
                n_bad_jobs += 1
            else:
                clear_logs(log_file.split(".")[0])

    print(f"Out of {n_total_jobs} jobs, {no_exit_code} had no exit code and {n_bad_jobs} had a bad exit code") 
    for code, info in exit_codes.items():
        print(f"Exit code {code} occurred {info['count']} times. Check {info['file']} for more info.")