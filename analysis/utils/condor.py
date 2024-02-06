import os

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
  err_files = [f for f in os.listdir(f"{condor_dir}/err") if job_tag in f]

  for err_file in err_files:
    with open(f"{condor_dir}/err/{err_file}") as f:
      # If the file is not empty
      if os.stat(f"{condor_dir}/err/{err_file}").st_size != 0:
        print(f"Error file {err_file} is not empty")
        print(err_file)